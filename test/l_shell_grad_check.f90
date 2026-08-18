!
! L shells through the gradient kernel: the packed answer against the split
! one, for int2e_ip1.
!
! Same contract as l_shell_check, which documents it in full: libcint cannot
! be given this input, so the reference is libfint's own split-shell path,
! which the gradient catalogue tests hold to the last bit against the C.
! The basis, the shared env, and the two error measures are l_shell_check's
! unchanged; what differs is the kernel under test and that every buffer
! carries the three derivative components.
!
! WHY THIS EXISTS AS A SEPARATE TEST.  A derivative integral is not just
! int2e with a different gout: it is the one case where the contraction
! buffers hold more than nf entries per primitive quartet.  The tensor
! component rides FASTEST there -- gout(n*3+m), transposed to
! component-slowest only at the very end -- so the stride at which an L
! shell's s/p coefficient pattern cycles carries a factor n_comp that a
! plain integral (n_comp = 1) never exercises.  That factor was missing
! when this test was first run: every third entry got the right
! coefficient and the rest were scrambled across components, half the
! tensor wrong by factors of hundreds.  l_shell_check passed throughout.
!
program l_shell_grad_check
   use cint_const, only: dp
   use cint_bas, only: ATM_SLOTS, BAS_SLOTS, ATOM_OF, ANG_OF, NPRIM_OF, &
                       NCTR_OF, KAPPA_OF, PTR_EXP, PTR_COEFF, &
                       CHARGE_OF, PTR_COORD, KAPPA_SP_SHELL, &
                       cint_gto_norm, cint_cgto_spheric, cint_cgto_cart
   use cint_envs, only: PTR_ENV_START
   use cint_workspace, only: cint_ws
   use cint_gen_grad2, only: int2e_ip1_sph, int2e_ip1_cart
   implicit none

   integer, parameter :: MAXBAS = 64, MAXENV = 4096, NATM = 3
   ! the three Cartesian components of nabla on centre i
   integer, parameter :: NCOMP = 3
   integer  :: atm(0:ATM_SLOTS*NATM-1)
   integer  :: pbas(0:BAS_SLOTS*MAXBAS-1), sbas(0:BAS_SLOTS*MAXBAS-1)
   real(dp) :: env(0:MAXENV-1)
   real(dp) :: coords(3, 0:NATM-1)
   integer  :: npbas, nsbas, envoff

   real(dp), allocatable :: eri_packed(:), eri_split(:)
   integer  :: nao_p, nao_s, sph
   real(dp) :: worst_rel, worst_abs
   integer  :: nbad
   real(dp), parameter :: TOL = 1.0e-13_dp
   real(dp), parameter :: TOL_ELEM = 1.0e-11_dp

   call build_basis()

   nbad = 0
   do sph = 1, 0, -1
      nao_p = nao_of(pbas, npbas, sph)
      nao_s = nao_of(sbas, nsbas, sph)
      if (nao_p /= nao_s) then
         print "(A,I0,A,I0)", "  FAIL: packed nao ", nao_p, " /= split nao ", nao_s
         stop 1
      end if
      allocate(eri_packed(0:NCOMP*nao_p**4-1), eri_split(0:NCOMP*nao_s**4-1))
      eri_packed = 0.0_dp; eri_split = 0.0_dp
      call all_eri(pbas, npbas, sph, nao_p, eri_packed)
      call all_eri(sbas, nsbas, sph, nao_s, eri_split)
      call compare(eri_packed, eri_split, NCOMP*nao_p**4, sph)
      deallocate(eri_packed, eri_split)
   end do

   if (nbad > 0) stop 1
   print "(A)", "  RESULT: PASS (int2e_ip1 L-packed reproduces split-shell)"

contains

   ! l_shell_check's basis, verbatim: an L shell of every shape the packing
   ! can produce -- segmented, generally contracted, one primitive -- plus a
   ! d shell so the cart-to-spherical stages run on quartets that mix an L
   ! shell with one that does transform.
   subroutine build_basis()
      real(dp) :: e3(3), c3(3), e2(2), c22(2,4), e1(1), c1(1)
      real(dp) :: cL3(3,2), cL1(1,2)

      atm = 0
      call set_atom(0, 6, [0.0_dp, 0.0_dp, 0.0_dp])
      call set_atom(1, 1, [0.0_dp, 0.0_dp, 1.9_dp])
      call set_atom(2, 8, [1.6_dp, 0.0_dp, -1.1_dp])

      npbas = 0; nsbas = 0
      envoff = PTR_ENV_START + 3*NATM
      env = 0.0_dp
      call set_atom_coords()

      ! carbon: a core s, a segmented L shell, a one-primitive L shell
      e3 = [71.6168370_dp, 13.0450960_dp, 3.5305122_dp]
      c3 = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
      call add_shell(0, 0, 3, 1, e3, reshape(c3, [3,1]), .false.)

      e3 = [2.9412494_dp, 0.6834831_dp, 0.2222899_dp]
      cL3(:,1) = [-0.09996723_dp, 0.39951283_dp, 0.70011547_dp]
      cL3(:,2) = [ 0.15591627_dp, 0.60768372_dp, 0.39195739_dp]
      call add_shell(0, 1, 3, 1, e3, cL3, .true.)

      e1 = [0.1687144_dp]
      cL1(:,1) = [1.0_dp]
      cL1(:,2) = [1.0_dp]
      call add_shell(0, 1, 1, 1, e1, cL1, .true.)

      ! hydrogen: two plain s shells
      e3 = [18.7311370_dp, 2.8253937_dp, 0.6401217_dp]
      c3 = [0.03349460_dp, 0.23472695_dp, 0.81375733_dp]
      call add_shell(1, 0, 3, 1, e3, reshape(c3, [3,1]), .false.)
      e1 = [0.1612778_dp]
      c1 = [1.0_dp]
      call add_shell(1, 0, 1, 1, e1, reshape(c1, [1,1]), .false.)

      ! oxygen: a core s, a GENERALLY CONTRACTED L shell -- two contractions
      ! sharing the same two primitives -- and a d polarisation shell
      e3 = [130.7093200_dp, 23.8088610_dp, 6.4436083_dp]
      c3 = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
      call add_shell(2, 0, 3, 1, e3, reshape(c3, [3,1]), .false.)

      e2 = [5.0331513_dp, 1.1695961_dp]
      c22(:,1) = [-0.09996723_dp, 0.39951283_dp]   ! s, contraction 0
      c22(:,2) = [ 0.21000000_dp, 0.83000000_dp]   ! s, contraction 1
      c22(:,3) = [ 0.15591627_dp, 0.60768372_dp]   ! p, contraction 0
      c22(:,4) = [ 0.44000000_dp,-0.27000000_dp]   ! p, contraction 1
      call add_shell(2, 1, 2, 2, e2, c22, .true.)

      e1 = [0.8000000_dp]
      c1 = [1.0_dp]
      call add_shell(2, 2, 1, 1, e1, reshape(c1, [1,1]), .false.)
   end subroutine build_basis

   subroutine set_atom(ia, z, r)
      integer,  intent(in) :: ia, z
      real(dp), intent(in) :: r(3)
      atm(ATM_SLOTS*ia + CHARGE_OF) = z
      atm(ATM_SLOTS*ia + PTR_COORD) = PTR_ENV_START + 3*ia
      coords(:, ia) = r
   end subroutine set_atom

   subroutine set_atom_coords()
      integer :: ia
      do ia = 0, NATM - 1
         env(PTR_ENV_START + 3*ia : PTR_ENV_START + 3*ia + 2) = coords(:, ia)
      end do
   end subroutine set_atom_coords

   ! One shell into the packed basis, and its equivalent into the split one.
   ! The split emission is s(ctr 0), p(ctr 0), s(ctr 1), p(ctr 1)... so the
   ! two bases share an AO ordering -- see l_shell_check for the argument.
   subroutine add_shell(ia, l, nprim, nctr, e, c, is_sp)
      integer,  intent(in) :: ia, l, nprim, nctr
      real(dp), intent(in) :: e(:), c(:,:)
      logical,  intent(in) :: is_sp
      integer :: ncol, m, ip, eptr, cptr, ic, lm

      ncol = nctr
      if (is_sp) ncol = 2*nctr
      eptr = envoff
      env(eptr:eptr+nprim-1) = e(1:nprim)
      envoff = envoff + nprim
      cptr = envoff
      do m = 1, ncol
         lm = l
         if (is_sp) then
            lm = 1
            if (m <= nctr) lm = 0
         end if
         do ip = 1, nprim
            env(cptr + (m-1)*nprim + ip-1) = c(ip, m) * cint_gto_norm(lm, e(ip))
         end do
      end do
      envoff = envoff + ncol*nprim

      call put(pbas, npbas, ia, l, nprim, nctr, eptr, cptr, is_sp)
      if (is_sp) then
         do ic = 0, nctr - 1
            call put(sbas, nsbas, ia, 0, nprim, 1, eptr, cptr + ic*nprim, .false.)
            call put(sbas, nsbas, ia, 1, nprim, 1, eptr, cptr + (nctr+ic)*nprim, .false.)
         end do
      else
         call put(sbas, nsbas, ia, l, nprim, nctr, eptr, cptr, .false.)
      end if
   end subroutine add_shell

   subroutine put(bas, nbas, ia, l, nprim, nctr, eptr, cptr, is_sp)
      integer, intent(inout) :: bas(0:), nbas
      integer, intent(in)    :: ia, l, nprim, nctr, eptr, cptr
      logical, intent(in)    :: is_sp
      integer :: b
      b = BAS_SLOTS*nbas
      bas(b + ATOM_OF)   = ia
      bas(b + ANG_OF)    = l
      bas(b + NPRIM_OF)  = nprim
      bas(b + NCTR_OF)   = nctr
      bas(b + KAPPA_OF)  = 0
      if (is_sp) bas(b + KAPPA_OF) = KAPPA_SP_SHELL
      bas(b + PTR_EXP)   = eptr
      bas(b + PTR_COEFF) = cptr
      nbas = nbas + 1
   end subroutine put

   integer function shell_dim(bas, ish, sph) result(d)
      integer, intent(in) :: bas(0:), ish, sph
      if (sph == 1) then
         d = cint_cgto_spheric(ish, bas)
      else
         d = cint_cgto_cart(ish, bas)
      end if
   end function shell_dim

   integer function nao_of(bas, nbas, sph) result(n)
      integer, intent(in) :: bas(0:), nbas, sph
      integer :: i
      n = 0
      do i = 0, nbas - 1
         n = n + shell_dim(bas, i, sph)
      end do
   end function nao_of

   ! Every quartet, scattered into one 3*nao^4 block -- the derivative
   ! component is the slowest index of both the buffer the kernel fills and
   ! the tensor compared, so the scatter carries it through unchanged.
   subroutine all_eri(bas, nbas, sph, nao, out)
      integer,  intent(in)    :: nbas, sph, nao
      integer,  intent(inout), target :: bas(0:)
      real(dp), intent(inout) :: out(0:)
      integer  :: ao_loc(0:MAXBAS), i, j, k, l, di, dj, dk, dl
      integer  :: shls(0:3), dims(0:3), a, b, c, d, m
      real(dp), allocatable :: buf(:)
      logical  :: hv
      type(cint_ws) :: ws

      ao_loc(0) = 0
      do i = 0, nbas - 1
         ao_loc(i+1) = ao_loc(i) + shell_dim(bas, i, sph)
      end do

      do i = 0, nbas - 1
      do j = 0, nbas - 1
      do k = 0, nbas - 1
      do l = 0, nbas - 1
         di = shell_dim(bas, i, sph); dj = shell_dim(bas, j, sph)
         dk = shell_dim(bas, k, sph); dl = shell_dim(bas, l, sph)
         allocate(buf(0:NCOMP*di*dj*dk*dl-1))
         buf = 0.0_dp
         shls = [i, j, k, l]
         dims = [di, dj, dk, dl]
         if (sph == 1) then
            hv = int2e_ip1_sph(buf, dims, shls, atm, NATM, bas, nbas, env, ws)
         else
            hv = int2e_ip1_cart(buf, dims, shls, atm, NATM, bas, nbas, env, ws)
         end if
         if (hv) then
            do m = 0, NCOMP - 1
            do d = 0, dl - 1
            do c = 0, dk - 1
            do b = 0, dj - 1
            do a = 0, di - 1
               out(m*nao**4 + ((ao_loc(l)+d)*nao + ao_loc(k)+c)*nao*nao &
                   + (ao_loc(j)+b)*nao + ao_loc(i)+a) &
                  = buf(m*di*dj*dk*dl + ((d*dk + c)*dj + b)*di + a)
            end do
            end do
            end do
            end do
            end do
         end if
         deallocate(buf)
      end do
      end do
      end do
      end do
   end subroutine all_eri

   ! Two measures, both reported -- l_shell_check documents why.  Scaled to
   ! the largest derivative in the set is the one the tolerance is on; per
   ! element over the entries above 1e-3 is what would catch a coefficient
   ! on the wrong sub-block or a component permuted, which is exactly the
   ! failure mode the n_comp stride bug produced.
   subroutine compare(a, b, n, sph)
      real(dp), intent(in) :: a(0:), b(0:)
      integer,  intent(in) :: n, sph
      integer  :: m, bad, nelem
      real(dp) :: d, r, scal, worst_elem
      character(len=4) :: what
      worst_rel = 0.0_dp; worst_abs = 0.0_dp; worst_elem = 0.0_dp
      bad = 0; nelem = 0
      scal = maxval(abs(b(0:n-1)))
      do m = 0, n - 1
         d = abs(a(m) - b(m))
         if (d > worst_abs) worst_abs = d
         r = d / scal
         if (r > worst_rel) worst_rel = r
         if (r > TOL) bad = bad + 1
         if (abs(b(m)) > 1.0e-3_dp) then
            nelem = nelem + 1
            r = d / abs(b(m))
            if (r > worst_elem) worst_elem = r
            if (r > TOL_ELEM) bad = bad + 1
         end if
      end do
      what = "cart"
      if (sph == 1) what = "sph "
      print "(A,A,1X,I0,A,ES10.2,A,ES10.2,A,I0,A,ES10.2,A,I0)", "  ", what, n, &
         "  worst abs ", worst_abs, "  scaled rel ", worst_rel, &
         "   |x|>1e-3: ", nelem, " worst rel ", worst_elem, "  over tol: ", bad
      nbad = nbad + bad
   end subroutine compare

end program l_shell_grad_check
