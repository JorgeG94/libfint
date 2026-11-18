!
! First derivatives of the scalar ECP integrals.
!
! Translated from PySCF's pyscf/lib/gto/nr_ecp_deriv.c -- the first-derivative
! part of it: `_deriv1_cart` and the four entry points built on it.
!
! What is differentiated is the bra.  A Cartesian Gaussian obeys
!
!     d/dA_x  x^l exp(-a r^2)  =  2a x^(l+1) exp(-a r^2)  -  l x^(l-1) exp(-a r^2)
!
! so the derivative of a shell of angular momentum l is a shell of l+1 and a
! shell of l-1, per primitive, and no new radial or angular code is needed:
! the kernels in cint_ecp are called twice on a fake, uncontracted shell pair
! and the two results are folded back.  The integrals returned are
!
!     ipnuc    < nabla_r mu | U | nu >          U summed over every ECP centre
!     iprinv   < nabla_r mu | U_C | nu >        U of one atom C alone
!
! three components each, x then y then z.  `nabla_r` is the derivative with
! respect to the electron's coordinate, which is minus the derivative with
! respect to the centre of `mu` -- libcint's `ip` convention, and the sign in
! the first line of `l_down`.
!
! The derivative with respect to the ECP centre is not here and is not
! needed: the integral is invariant to translating all three centres
! together, so it is minus the sum of the two basis-centre derivatives, and
! `iprinv` is what restricts that sum to one atom's potential.
!
! The atom `iprinv` refers to is read from `env(AS_RINV_ORIG_ATOM)`, zero
! based, which is the slot PySCF's `with_rinv_at_nucleus` writes.  The ECP
! shells of one atom have to be contiguous in `ecpbas`, which is how PySCF
! lays them out.
!
! The Hessian's integrals (ipipnuc, ipnucip and their one-atom forms) are in
! the same C file and are not translated.
!
module cint_ecp_deriv
   use cint_const, only: dp
   use cint_bas, only: ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, &
                       BAS_SLOTS, cint_len_cart
   use cint_ecp_num, only: AS_RINV_ORIG_ATOM, AS_ECPBAS_OFFSET, AS_NECPBAS
   use cint_ecp, only: ecp_type1_cart, ecp_type2_cart
   use cint_ecp_drv, only: ecp_cart_to_sph, ecp_distribute
   implicit none
   private

   public :: ecp_ipnuc_cart, ecp_ipnuc_sph
   public :: ecp_iprinv_cart, ecp_iprinv_sph
   public :: ecp_deriv1, ecp_deriv1_cart

   !> Components of a first derivative
   integer, parameter :: NCOMP = 3

   !> Highest total angular momentum the angular code's tables hold: CART_CUM
   !> in cint_ecp_ang is the count of Cartesian functions through l = 12.
   integer, parameter :: ANG_LMAX = 12

   !> The factor the kernels fold into a p shell, sqrt(3/(4 pi))
   real(dp), parameter :: FAC_P = 0.488602511902919921_dp

contains

   !> Where one atom's ECP shells sit in `ecpbas`
   !>
   !> `first` is the zero-based row of the first shell on `atm_id`, or -1 when
   !> that atom carries no ECP; `nsh` is how many consecutive rows belong to it.
   pure subroutine one_atom_ecpbas(atm_id, ecpbas, necpbas, first, nsh)
      integer, intent(in) :: atm_id, ecpbas(0:), necpbas
      integer, intent(out) :: first, nsh
      integer :: i

      first = -1
      nsh = 0
      do i = 0, necpbas - 1
         if (ecpbas(ATOM_OF + i*BAS_SLOTS) == atm_id) then
            if (first < 0) first = i
            nsh = nsh + 1
         else if (first >= 0) then
            exit
         end if
      end do
   end subroutine one_atom_ecpbas

   !> The shell pair as primitives: `npi` rows for the bra, then `npj` for the ket
   !>
   !> Each row is one uncontracted primitive whose "coefficient" pointer is
   !> its own exponent, because a pointer into `env` is all a `bas` row can
   !> carry and the exponent is a number known to be there.  The caller
   !> divides it back out.
   pure subroutine uncontract_bas(fakbas, shls, bas)
      integer, intent(out) :: fakbas(0:)
      integer, intent(in) :: shls(0:1), bas(0:)
      integer :: side, sh, np, i, row

      fakbas = 0
      row = 0
      do side = 0, 1
         sh = shls(side)
         np = bas(NPRIM_OF + sh*BAS_SLOTS)
         do i = 0, np - 1
            fakbas(ATOM_OF + row*BAS_SLOTS) = bas(ATOM_OF + sh*BAS_SLOTS)
            fakbas(ANG_OF + row*BAS_SLOTS) = bas(ANG_OF + sh*BAS_SLOTS)
            fakbas(NPRIM_OF + row*BAS_SLOTS) = 1
            fakbas(NCTR_OF + row*BAS_SLOTS) = 1
            fakbas(PTR_EXP + row*BAS_SLOTS) = bas(PTR_EXP + sh*BAS_SLOTS) + i
            fakbas(PTR_COEFF + row*BAS_SLOTS) = bas(PTR_EXP + sh*BAS_SLOTS) + i
            row = row + 1
         end do
      end do
   end subroutine uncontract_bas

   !> The `-2a x^(l+1)` half: read a derivative off the raised shell
   !>
   !> `buf1` is `(nfi1, nfj)` over the shell of `li + 1`; `out` is three
   !> `(nfi, nfj)` blocks and is overwritten.  A Cartesian function at
   !> position `i` of a shell, with `m = ly + lz`, has its x-raised partner at
   !> the same position of the next shell, its y-raised one at `i + m + 1`
   !> and its z-raised one at `i + m + 2`.
   !>
   !> The kernels fold a spherical-harmonic factor into s and p shells and
   !> none into the rest, so a derivative that crosses from s to p or from p
   !> to d has to take one factor out and put the other in.
   pure subroutine l_down(out, buf1, fac, ai, li, nfj)
      real(dp), intent(out) :: out(0:)
      real(dp), intent(in) :: buf1(0:), fac, ai
      integer, intent(in) :: li, nfj

      integer :: nfi, nfi1, i, j, m, lz, nij
      real(dp) :: f

      nfi = cint_len_cart(li)
      nfi1 = cint_len_cart(li + 1)
      nij = nfi*nfj

      select case (li)
      case (0)
         f = fac*(-2.0_dp/sqrt(3.0_dp))*ai
      case (1)
         f = fac*(-2.0_dp*FAC_P)*ai
      case default
         f = fac*(-2.0_dp)*ai
      end select

      i = 0
      do m = 0, li
         do lz = 0, m
            do j = 0, nfj - 1
               out(j*nfi + i) = f*buf1(j*nfi1 + i)
               out(nij + j*nfi + i) = f*buf1(j*nfi1 + i + m + 1)
               out(2*nij + j*nfi + i) = f*buf1(j*nfi1 + i + m + 2)
            end do
            i = i + 1
         end do
      end do
   end subroutine l_down

   !> The `l x^(l-1)` half: add a derivative read off the lowered shell
   !>
   !> `buf1` is `(nfi0, nfj)` over the shell of `li - 1`; `out` is accumulated
   !> into.  Function `i` of the lowered shell, with powers `(lx, ly, lz)`,
   !> contributes to the x derivative of its x-raised partner with weight
   !> `lx + 1`, and likewise for y and z.
   pure subroutine l_up(out, buf1, fac, li, nfj)
      real(dp), intent(inout) :: out(0:)
      real(dp), intent(in) :: buf1(0:), fac
      integer, intent(in) :: li, nfj

      integer :: nfi, nfi0, i, j, m, lz, ly, lx, nij
      real(dp) :: f

      nfi = cint_len_cart(li)
      nfi0 = cint_len_cart(li - 1)
      nij = nfi*nfj

      select case (li)
      case (1)
         f = fac*sqrt(3.0_dp)
      case (2)
         f = fac/FAC_P
      case default
         f = fac
      end select

      i = 0
      do m = 0, li - 1
         do lz = 0, m
            ly = m - lz
            lx = li - 1 - m
            do j = 0, nfj - 1
               out(j*nfi + i) = out(j*nfi + i) + f*(lx + 1)*buf1(j*nfi0 + i)
               out(nij + j*nfi + i + m + 1) = out(nij + j*nfi + i + m + 1) &
                                              + f*(ly + 1)*buf1(j*nfi0 + i)
               out(2*nij + j*nfi + i + m + 2) = out(2*nij + j*nfi + i + m + 2) &
                                                + f*(lz + 1)*buf1(j*nfi0 + i)
            end do
            i = i + 1
         end do
      end do
   end subroutine l_up

   !> `< nabla mu | U | nu >` over one Cartesian shell pair, three components
   !>
   !> Accumulates into `gctr`, three `(nfi*nci, nfj*ncj)` blocks, which the
   !> caller has zeroed.  `ecpbas` is whichever shells make up `U`: all of
   !> them, or one atom's.  Returns non-zero when anything was added.
   function ecp_deriv1_cart(gctr, shls, ecpbas, necpbas, atm, natm, bas, env) &
      result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer, intent(in) :: shls(0:1), ecpbas(0:), necpbas, atm(0:), bas(0:), natm
      real(dp), intent(in) :: env(0:)
      integer :: has_value

      integer :: ish, jsh, li, lj, npi, npj, nci, ncj, nfi, nfj, nfi1, di, dj, dij
      integer :: pai, paj, pci, pcj, nfak, ip, jp, ic, jc, i, j, n, c, shls1(0:1)
      real(dp) :: fac
      integer, allocatable :: fakbas(:)
      real(dp), allocatable :: buf1(:), gprim(:)

      has_value = 0
      if (necpbas == 0) return

      ish = shls(0); jsh = shls(1)
      li = bas(ANG_OF + ish*BAS_SLOTS); lj = bas(ANG_OF + jsh*BAS_SLOTS)
      ! The raised bra takes the angular code to li + lj + 1, and its tables
      ! stop at ANG_LMAX.  Only a pair of i shells gets there; the C has no
      ! guard and writes past its tables for one.
      if (li + lj + 1 > ANG_LMAX) then
         error stop "cint_ecp_deriv: li + lj + 1 exceeds the ECP angular tables (two i shells)"
      end if
      npi = bas(NPRIM_OF + ish*BAS_SLOTS); npj = bas(NPRIM_OF + jsh*BAS_SLOTS)
      nci = bas(NCTR_OF + ish*BAS_SLOTS); ncj = bas(NCTR_OF + jsh*BAS_SLOTS)
      nfi = cint_len_cart(li); nfj = cint_len_cart(lj)
      nfi1 = cint_len_cart(li + 1)
      di = nfi*nci; dj = nfj*ncj; dij = di*dj
      pai = bas(PTR_EXP + ish*BAS_SLOTS); paj = bas(PTR_EXP + jsh*BAS_SLOTS)
      pci = bas(PTR_COEFF + ish*BAS_SLOTS); pcj = bas(PTR_COEFF + jsh*BAS_SLOTS)

      nfak = npi + npj
      allocate (fakbas(0:nfak*BAS_SLOTS - 1))
      call uncontract_bas(fakbas, shls, bas)
      allocate (buf1(0:nfi1*nfj - 1), gprim(0:NCOMP*nfi*nfj - 1))

      do jp = 0, npj - 1
         do ip = 0, npi - 1
            shls1 = [ip, npi + jp]
            ! The fake shells carry their exponents as coefficients.
            fac = 1.0_dp/(env(pai + ip)*env(paj + jp))

            fakbas(ANG_OF + ip*BAS_SLOTS) = li + 1
            buf1 = 0.0_dp
            has_value = ior(has_value, ecp_type1_cart(buf1, shls1, ecpbas, necpbas, &
                                                      atm, natm, fakbas, nfak, env))
            has_value = ior(has_value, ecp_type2_cart(buf1, shls1, ecpbas, necpbas, &
                                                      atm, natm, fakbas, nfak, env))
            call l_down(gprim, buf1, fac, env(pai + ip), li, nfj)

            if (li > 0) then
               fakbas(ANG_OF + ip*BAS_SLOTS) = li - 1
               buf1 = 0.0_dp
               has_value = ior(has_value, ecp_type1_cart(buf1, shls1, ecpbas, necpbas, &
                                                         atm, natm, fakbas, nfak, env))
               has_value = ior(has_value, ecp_type2_cart(buf1, shls1, ecpbas, necpbas, &
                                                         atm, natm, fakbas, nfak, env))
               call l_up(gprim, buf1, fac, li, nfj)
            end if
            fakbas(ANG_OF + ip*BAS_SLOTS) = li

            do jc = 0, ncj - 1
               do ic = 0, nci - 1
                  fac = env(pci + ic*npi + ip)*env(pcj + jc*npj + jp)
                  n = jc*nfj*di + ic*nfi
                  do c = 0, NCOMP - 1
                     do j = 0, nfj - 1
                        do i = 0, nfi - 1
                           gctr(c*dij + n + j*di + i) = gctr(c*dij + n + j*di + i) &
                                                        + fac*gprim(c*nfi*nfj + j*nfi + i)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end function ecp_deriv1_cart

   !> The shell pair's block in the requested basis, placed in `out`
   !>
   !> The potential is `nsh` rows of `ecpbas` starting at row `first`.  With
   !> `first < 0` or `nsh == 0` it is empty: `out` is zeroed and the result is
   !> 0, as for a pair the overlap screen rejects.
   function deriv1_driver(cart, out, shls, first, nsh, atm, natm, bas, env, dims) &
      result(has_value)
      logical, intent(in) :: cart
      real(dp), intent(inout) :: out(0:)
      integer, intent(in) :: shls(0:1), first, nsh, atm(0:), bas(0:), natm
      real(dp), intent(in) :: env(0:)
      integer, intent(in), optional :: dims(0:1)
      integer :: has_value

      integer :: li, lj, nci, ncj, dic, djc, di, dj, c, nout, offset
      real(dp), allocatable :: gcart(:), gsph(:)

      li = bas(ANG_OF + shls(0)*BAS_SLOTS); lj = bas(ANG_OF + shls(1)*BAS_SLOTS)
      nci = bas(NCTR_OF + shls(0)*BAS_SLOTS); ncj = bas(NCTR_OF + shls(1)*BAS_SLOTS)
      dic = cint_len_cart(li)*nci; djc = cint_len_cart(lj)*ncj
      if (cart) then
         di = dic; dj = djc
      else
         di = (2*li + 1)*nci; dj = (2*lj + 1)*ncj
      end if
      nout = di*dj
      if (present(dims)) nout = dims(0)*dims(1)

      allocate (gcart(0:NCOMP*dic*djc - 1))
      gcart = 0.0_dp
      has_value = 0
      if (first >= 0 .and. nsh > 0) then
         offset = int(env(AS_ECPBAS_OFFSET))*BAS_SLOTS + first*BAS_SLOTS
         has_value = ecp_deriv1_cart(gcart, shls, bas(offset:), nsh, atm, natm, bas, env)
      end if

      if (cart) then
         do c = 0, NCOMP - 1
            call ecp_distribute(out(c*nout:), gcart(c*dic*djc:), dims, di, dj, has_value == 0)
         end do
      else
         allocate (gsph(0:di*dj - 1))
         do c = 0, NCOMP - 1
            if (has_value /= 0) call ecp_cart_to_sph(gsph, gcart(c*dic*djc:), li, lj, nci, ncj)
            call ecp_distribute(out(c*nout:), gsph, dims, di, dj, has_value == 0)
         end do
      end if
   end function deriv1_driver

   !> The ECP rows of the atom `env` names
   pure subroutine atom_shells(env, bas, first, nsh)
      real(dp), intent(in) :: env(0:)
      integer, intent(in) :: bas(0:)
      integer, intent(out) :: first, nsh
      integer :: offset, necpbas

      first = -1
      nsh = 0
      necpbas = int(env(AS_NECPBAS))
      if (necpbas <= 0) return
      offset = int(env(AS_ECPBAS_OFFSET))*BAS_SLOTS
      call one_atom_ecpbas(int(env(AS_RINV_ORIG_ATOM)), bas(offset:), necpbas, first, nsh)
   end subroutine atom_shells

   !> Either first-derivative integral over a shell pair, in either basis
   !>
   !> `one_atom` selects `iprinv`, the potential of the atom in
   !> `env(AS_RINV_ORIG_ATOM)` alone, over `ipnuc`, every ECP centre.  `out`
   !> holds three blocks, x then y then z, each `(di, dj)` -- or each
   !> `dims(0)*dims(1)` with the block in its leading corner when `dims` is
   !> present.  The result is 0, and `out` is zeroed, when nothing
   !> contributes: no ECP, an atom without one, or a pair the screen rejects.
   function ecp_deriv1(cart, one_atom, out, shls, atm, natm, bas, env, dims) &
      result(has_value)
      logical, intent(in) :: cart, one_atom
      real(dp), intent(inout) :: out(0:)
      integer, intent(in) :: shls(0:1), atm(0:), bas(0:), natm
      real(dp), intent(in) :: env(0:)
      integer, intent(in), optional :: dims(0:1)
      integer :: has_value
      integer :: first, nsh

      if (one_atom) then
         call atom_shells(env, bas, first, nsh)
      else
         first = 0
         nsh = int(env(AS_NECPBAS))
      end if
      has_value = deriv1_driver(cart, out, shls, first, nsh, atm, natm, bas, env, dims)
   end function ecp_deriv1

   !> `< nabla mu | U | nu >` over a shell pair, spherical basis
   !>
   !> See `ecp_deriv1`.  `nbas` is not read: the ECP rows are found through
   !> `env`.  It is in the signature because every other integral here takes it.
   function ecp_ipnuc_sph(out, shls, atm, natm, bas, nbas, env, dims) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer, intent(in) :: shls(0:1), atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer, intent(in), optional :: dims(0:1)
      integer :: has_value

      has_value = ecp_deriv1(.false., .false., out, shls, atm, natm, bas, env, dims)
   end function ecp_ipnuc_sph

   !> `< nabla mu | U | nu >` over a shell pair, Cartesian basis
   !>
   !> See `ecp_deriv1`.  `nbas` is not read: the ECP rows are found through
   !> `env`.  It is in the signature because every other integral here takes it.
   function ecp_ipnuc_cart(out, shls, atm, natm, bas, nbas, env, dims) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer, intent(in) :: shls(0:1), atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer, intent(in), optional :: dims(0:1)
      integer :: has_value

      has_value = ecp_deriv1(.true., .false., out, shls, atm, natm, bas, env, dims)
   end function ecp_ipnuc_cart

   !> `< nabla mu | U_C | nu >` for the atom in `env(AS_RINV_ORIG_ATOM)`, spherical
   !>
   !> See `ecp_deriv1`.  `nbas` is not read: the ECP rows are found through
   !> `env`.  It is in the signature because every other integral here takes it.
   function ecp_iprinv_sph(out, shls, atm, natm, bas, nbas, env, dims) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer, intent(in) :: shls(0:1), atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer, intent(in), optional :: dims(0:1)
      integer :: has_value

      has_value = ecp_deriv1(.false., .true., out, shls, atm, natm, bas, env, dims)
   end function ecp_iprinv_sph

   !> `< nabla mu | U_C | nu >` for the atom in `env(AS_RINV_ORIG_ATOM)`, Cartesian
   !>
   !> See `ecp_deriv1`.  `nbas` is not read: the ECP rows are found through
   !> `env`.  It is in the signature because every other integral here takes it.
   function ecp_iprinv_cart(out, shls, atm, natm, bas, nbas, env, dims) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer, intent(in) :: shls(0:1), atm(0:), bas(0:), natm, nbas
      real(dp), intent(in) :: env(0:)
      integer, intent(in), optional :: dims(0:1)
      integer :: has_value

      has_value = ecp_deriv1(.true., .true., out, shls, atm, natm, bas, env, dims)
   end function ecp_iprinv_cart

end module cint_ecp_deriv
