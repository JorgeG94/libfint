!
! Basis-set bookkeeping: how many components a shell has, where each shell
! starts in the AO list, and the Cartesian component ordering.
!
! Ported from src/cint_bas.c and the numeric helpers of src/misc.c.
!
! INDEXING.  `bas` and `atm` keep the C's flat layout and 0-based addressing:
! the C spells it with a macro, `bas(SLOT,I) -> bas[BAS_SLOTS*I + SLOT]`, and
! that arithmetic is reproduced here rather than reshaped into a 2D array.  The
! slot constants are therefore the C's 0-based values, NOT the 1-based ones in
! libcint_interface.f90 -- the two conventions differ and confusing them is
! precisely the bug D0 found in the timing benchmark.
!
module cint_bas
   use cint_const, only: dp
   implicit none
   private

   public :: cint_len_cart, cint_len_spinor
   public :: cint_cgto_cart, cint_cgto_spheric, cint_cgto_spinor
   public :: cint_tot_pgto_spheric, cint_tot_pgto_spinor
   public :: cint_tot_cgto_cart, cint_tot_cgto_spheric, cint_tot_cgto_spinor
   public :: cint_shells_cart_offset, cint_shells_spheric_offset
   public :: cint_shells_spinor_offset
   public :: cint_cart_comp, cint_cart_comp_sp
   public :: cint_square_dist, cint_gto_norm
   public :: cint_bas_is_sp

   ! Slot offsets, 0-based, as in include/cint.h.in.
   integer, parameter, public :: ATOM_OF   = 0
   integer, parameter, public :: ANG_OF    = 1
   integer, parameter, public :: NPRIM_OF  = 2
   integer, parameter, public :: NCTR_OF   = 3
   integer, parameter, public :: KAPPA_OF  = 4
   integer, parameter, public :: PTR_EXP   = 5
   integer, parameter, public :: PTR_COEFF = 6
   integer, parameter, public :: BAS_SLOTS = 8

   ! L SHELLS.  A shell that carries both s and p on one set of exponents,
   ! with a separate contraction column for each -- what GAMESS calls an L
   ! shell and what every Pople basis is written in.  `bas` has one ANG_OF
   ! slot, so libcint cannot say it and splits such a shell in two at the
   ! packing boundary; this marker says it instead.
   !
   ! The convention, and it is checked nowhere because `bas` is the caller's:
   !
   !   ANG_OF    = 1          -- the highest l present, so every ceiling,
   !                             stride and Rys order the machinery derives
   !                             is the one the p sub-block needs.  The s
   !                             sub-block is then computed with more roots
   !                             than it needs, which is exact: Rys is exact
   !                             to degree 2n-1.
   !   KAPPA_OF  = KAPPA_SP_SHELL
   !   PTR_COEFF = nprim*nctr s coefficients, then nprim*nctr p ones -- the
   !               same block a plain shell with 2*nctr contractions would
   !               have, which is why no packer has to move anything.
   !
   ! KAPPA_OF is the relativistic kappa, zero on every non-relativistic
   ! shell and within +-(ANG_MAX+1) on every spinor one, so a value outside
   ! that range cannot collide with a real kappa.
   integer, parameter, public :: KAPPA_SP_SHELL = 64
   ! s plus p: one component then three, in that order.
   integer, parameter, public :: NF_SP = 4

   integer, parameter, public :: CHARGE_OF       = 0
   integer, parameter, public :: PTR_COORD       = 1
   integer, parameter, public :: NUC_MOD_OF      = 2
   integer, parameter, public :: PTR_ZETA        = 3
   integer, parameter, public :: PTR_FRAC_CHARGE = 4
   integer, parameter, public :: ATM_SLOTS       = 6

   ! How a nucleus is modelled: a point, a Gaussian of width zeta, or a
   ! fractional charge.
   integer, parameter, public :: POINT_NUC       = 1
   integer, parameter, public :: GAUSSIAN_NUC    = 2
   integer, parameter, public :: FRAC_CHARGE_NUC = 3

contains

   ! bas(SLOT, I) in the C's spelling.
   pure function basval(bas, slot, i) result(v)
      integer, intent(in) :: bas(0:), slot, i
      integer :: v
      v = bas(BAS_SLOTS * i + slot)
   end function basval

   ! Components of a Cartesian GTO shell, (l+1)(l+2)/2.
   pure function cint_len_cart(l) result(n)
      integer, intent(in) :: l
      integer :: n
      n = (l + 1) * (l + 2) / 2
   end function cint_len_cart

   ! Components of a spinor shell, which depends on kappa: zero means both
   ! j = l +- 1/2 are present, negative means j = l + 1/2 alone, positive
   ! means j = l - 1/2 alone.
   pure function cint_len_spinor(bas_id, bas) result(n)
      integer, intent(in) :: bas_id, bas(0:)
      integer :: n, kappa, l
      kappa = basval(bas, KAPPA_OF, bas_id)
      l = basval(bas, ANG_OF, bas_id)
      if (kappa == 0) then
         n = 4 * l + 2
      else if (kappa < 0) then
         n = 2 * l + 2
      else
         n = 2 * l
      end if
   end function cint_len_spinor

   ! Does this shell carry both s and p on one set of exponents?  See the
   ! note on KAPPA_SP_SHELL above.
   pure function cint_bas_is_sp(bas_id, bas) result(yes)
      integer, intent(in) :: bas_id, bas(0:)
      logical :: yes
      yes = basval(bas, KAPPA_OF, bas_id) == KAPPA_SP_SHELL
   end function cint_bas_is_sp

   pure function cint_cgto_cart(bas_id, bas) result(n)
      integer, intent(in) :: bas_id, bas(0:)
      integer :: n, l
      if (cint_bas_is_sp(bas_id, bas)) then
         n = NF_SP * basval(bas, NCTR_OF, bas_id)
         return
      end if
      l = basval(bas, ANG_OF, bas_id)
      n = (l+1)*(l+2)/2 * basval(bas, NCTR_OF, bas_id)
   end function cint_cgto_cart

   ! s and p are the two l for which the spherical functions ARE the
   ! Cartesian ones, in the same order -- which is why cint_c2s_*_sph return
   ! RESULT_IN_GCART below l = 2.  An L shell is therefore four functions
   ! either way, and needs no transform at all.
   pure function cint_cgto_spheric(bas_id, bas) result(n)
      integer, intent(in) :: bas_id, bas(0:)
      integer :: n
      if (cint_bas_is_sp(bas_id, bas)) then
         n = NF_SP * basval(bas, NCTR_OF, bas_id)
         return
      end if
      n = (basval(bas, ANG_OF, bas_id) * 2 + 1) * basval(bas, NCTR_OF, bas_id)
   end function cint_cgto_spheric

   pure function cint_cgto_spinor(bas_id, bas) result(n)
      integer, intent(in) :: bas_id, bas(0:)
      integer :: n
      n = cint_len_spinor(bas_id, bas) * basval(bas, NCTR_OF, bas_id)
   end function cint_cgto_spinor

   pure function cint_tot_pgto_spheric(bas, nbas) result(s)
      integer, intent(in) :: bas(0:), nbas
      integer :: s, i
      s = 0
      do i = 0, nbas - 1
         s = s + (basval(bas, ANG_OF, i) * 2 + 1) * basval(bas, NPRIM_OF, i)
      end do
   end function cint_tot_pgto_spheric

   pure function cint_tot_pgto_spinor(bas, nbas) result(s)
      integer, intent(in) :: bas(0:), nbas
      integer :: s, i
      s = 0
      do i = 0, nbas - 1
         s = s + cint_len_spinor(i, bas) * basval(bas, NPRIM_OF, i)
      end do
   end function cint_tot_pgto_spinor

   ! The C reaches these three through a function pointer (tot_cgto_accum
   ! taking `FINT (*f)()`).  Three short loops say the same thing without an
   ! abstract interface for a one-line body.
   pure function cint_tot_cgto_cart(bas, nbas) result(s)
      integer, intent(in) :: bas(0:), nbas
      integer :: s, i
      s = 0
      do i = 0, nbas - 1
         s = s + cint_cgto_cart(i, bas)
      end do
   end function cint_tot_cgto_cart

   pure function cint_tot_cgto_spheric(bas, nbas) result(s)
      integer, intent(in) :: bas(0:), nbas
      integer :: s, i
      s = 0
      do i = 0, nbas - 1
         s = s + cint_cgto_spheric(i, bas)
      end do
   end function cint_tot_cgto_spheric

   pure function cint_tot_cgto_spinor(bas, nbas) result(s)
      integer, intent(in) :: bas(0:), nbas
      integer :: s, i
      s = 0
      do i = 0, nbas - 1
         s = s + cint_cgto_spinor(i, bas)
      end do
   end function cint_tot_cgto_spinor

   ! ao_loc(0:nbas) -- the offset of each shell in the AO list, plus the total.
   pure subroutine cint_shells_cart_offset(ao_loc, bas, nbas)
      integer, intent(out) :: ao_loc(0:)
      integer, intent(in)  :: bas(0:), nbas
      integer :: i
      ao_loc(0) = 0
      do i = 0, nbas - 1
         ao_loc(i+1) = ao_loc(i) + cint_cgto_cart(i, bas)
      end do
   end subroutine cint_shells_cart_offset

   pure subroutine cint_shells_spheric_offset(ao_loc, bas, nbas)
      integer, intent(out) :: ao_loc(0:)
      integer, intent(in)  :: bas(0:), nbas
      integer :: i
      ao_loc(0) = 0
      do i = 0, nbas - 1
         ao_loc(i+1) = ao_loc(i) + cint_cgto_spheric(i, bas)
      end do
   end subroutine cint_shells_spheric_offset

   pure subroutine cint_shells_spinor_offset(ao_loc, bas, nbas)
      integer, intent(out) :: ao_loc(0:)
      integer, intent(in)  :: bas(0:), nbas
      integer :: i
      ao_loc(0) = 0
      do i = 0, nbas - 1
         ao_loc(i+1) = ao_loc(i) + cint_cgto_spinor(i, bas)
      end do
   end subroutine cint_shells_spinor_offset

   ! Cartesian component powers for a shell of angular momentum lmax, in
   ! libcint's order: x descending outermost, then y.  This order is a
   ! convention the rest of the library depends on, so it is reproduced
   ! exactly rather than derived.
   pure subroutine cint_cart_comp(nx, ny, nz, lmax)
      integer, intent(out) :: nx(0:), ny(0:), nz(0:)
      integer, intent(in)  :: lmax
      integer :: inc, lx, ly, lz
      inc = 0
      do lx = lmax, 0, -1
         do ly = lmax - lx, 0, -1
            lz = lmax - lx - ly
            nx(inc) = lx
            ny(inc) = ly
            nz(inc) = lz
            inc = inc + 1
         end do
      end do
   end subroutine cint_cart_comp

   ! The same, for an L shell: the s function first, then the three p in the
   ! order a plain l = 1 shell has them.  That order is what makes the sp
   ! index map below the l = 1 map with one entry prepended, and what lets
   ! the cart-to-spherical stage stay the identity.
   pure subroutine cint_cart_comp_sp(nx, ny, nz)
      integer, intent(out) :: nx(0:), ny(0:), nz(0:)
      nx(0:3) = [0, 1, 0, 0]
      ny(0:3) = [0, 0, 1, 0]
      nz(0:3) = [0, 0, 0, 1]
   end subroutine cint_cart_comp_sp

   pure function cint_square_dist(r1, r2) result(d)
      real(dp), intent(in) :: r1(0:), r2(0:)
      real(dp) :: d, r12(0:2)
      r12(0) = r1(0) - r2(0)
      r12(1) = r1(1) - r2(1)
      r12(2) = r1(2) - r2(2)
      d = r12(0)*r12(0) + r12(1)*r12(1) + r12(2)*r12(2)
   end function cint_square_dist

   ! \int_0^\infty r^n e^{-alpha r^2} dr
   pure function gaussian_int(n, alpha) result(g)
      integer,  intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp) :: g, n1
      n1 = (n + 1) * 0.5_dp
      g = exp(log_gamma(n1)) / (2.0_dp * alpha**n1)
   end function gaussian_int

   ! Normalisation of the GTO radial part g = r^l e^{-a r^2}, i.e.
   ! 1/sqrt(int g^2 r^2 dr).
   ! Ref: Schlegel and Frisch, Int. J. Quant. Chem. 54 (1995) 83.
   pure function cint_gto_norm(n, a) result(f)
      integer,  intent(in) :: n
      real(dp), intent(in) :: a
      real(dp) :: f
      f = 1.0_dp / sqrt(gaussian_int(n*2 + 2, 2.0_dp*a))
   end function cint_gto_norm

end module cint_bas
