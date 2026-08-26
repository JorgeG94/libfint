!
! Angular factors for the ECP integrals.
!
! Translated from PySCF's pyscf/lib/gto/nr_ecp.c.  See cint_ecp_num for why
! the upstream here is PySCF rather than libcint.
!
! Both ECP integral types reduce to a radial quadrature multiplied by an
! angular factor that depends only on the shell's angular momenta and the
! direction from the ECP centre.  This module builds those factors.  Nothing
! here touches a primitive exponent, which is what lets the type-2 factors be
! computed once per shell pair and reused across every primitive.
!
! The one thing worth knowing before reading it: the projection through
! spherical harmonics in `ang_nuc_in_cart` is what makes the angular integral
! exact.  A Cartesian monomial of degree l is not a pure spherical harmonic --
! it carries lower-l contamination -- and the round trip through the
! spherical basis and back is what removes it.
!
module cint_ecp_ang
   use cint_const, only: dp
   use cint_bas, only: cint_len_cart
   use cint_cart2sph, only: cint_c2s_bra_sph, cint_s2c_bra_sph, RESULT_IN_GCART
   use cint_ecp_num, only: ecp_factorial2, ecp_binom
   implicit none
   private

   public :: ecp_ang_nuc_in_cart
   public :: ecp_int_unit_xyz
   public :: ecp_cache_3dfac
   public :: ecp_type2_facs_ang
   public :: ecp_type1_rad_ang
   public :: CART_POW_Y, CART_POW_Z, OFFSET_CART, CART_CUM

   !> Cumulative Cartesian component counts, so shell l starts at OFFSET_CART(l)
   integer, parameter :: OFFSET_CART(0:14) = &
                         [0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 364, 455, 560]

   !> Total Cartesian components up to l = 12, the size of an omega_nuc buffer
   integer, parameter :: CART_CUM = 456

   ! The y and z powers of the Cartesian components, in libcint's ordering.
   !
   ! Indexed from zero **regardless of l**: the C's LOOP_CART starts its
   ! counter at 0 and not at OFFSET_CART(l), so the first (l+1)(l+2)/2 entries
   ! of these tables are shell l's ordering. Reading them at OFFSET_CART(l)
   ! instead -- which is what the name suggests and what OFFSET_CART is for
   ! elsewhere in this file -- silently produces a different shell's exponents.
   integer, parameter :: CART_POW_Y(0:119) = [ &
                         0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, &
                         0, 6, 5, 4, 3, 2, 1, 0, 7, 6, 5, 4, 3, 2, 1, 0, 8, 7, 6, 5, &
                         4, 3, 2, 1, 0, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 9, 8, 7, 6, &
                         5, 4, 3, 2, 1, 0, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 12, 11, &
                         10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 13, 12, 11, 10, 9, 8, 7, 6, 5, &
                         4, 3, 2, 1, 0, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]

   integer, parameter :: CART_POW_Z(0:119) = [ &
                         0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, &
                         5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, &
                         4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, &
                         5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, &
                         2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8, &
                         9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

contains

   !> The unit vector from a displacement, or zero when there is no direction
   !>
   !> Negated, as the C has it: the angular factors are built about the vector
   !> from the shell centre towards the ECP, and that sign is carried here
   !> rather than at each of the several call sites.
   pure subroutine unit_direction(r, unitr)
      real(dp), intent(in) :: r(3)
      real(dp), intent(out) :: unitr(3)
      real(dp) :: r2, norm

      r2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
      if (r(1) == 0.0_dp .and. r(2) == 0.0_dp .and. r(3) == 0.0_dp) then
         unitr = 0.0_dp
      else
         norm = -1.0_dp/sqrt(r2)
         unitr = r*norm
      end if
   end subroutine unit_direction

   !> The pure solid harmonic of degree l, written in Cartesian components
   !>
   !> For l >= 2 a Cartesian monomial of degree l is not a pure harmonic, so
   !> the monomials are projected onto the spherical basis and transformed
   !> straight back. That round trip is the projection: what survives it is
   !> the l-pure part, which is the only part that contributes to the angular
   !> integral.
   !>
   !> l = 0 and l = 1 are closed forms because the transform is the identity
   !> there, and the constants are 1/(4 pi) and 3/(4 pi) folded in.
   subroutine ecp_ang_nuc_in_cart(omega, l, r)
      real(dp), intent(inout) :: omega(0:)
      integer, intent(in) :: l
      real(dp), intent(in) :: r(3)

      real(dp) :: buf(0:(14 + 1)*(14 + 2)/2 - 1)
      real(dp) :: xx(0:15), yy(0:15), zz(0:15)
      integer :: i, j, k, n, loc

      select case (l)
      case (0)
         omega(0) = 0.07957747154594767_dp
      case (1)
         omega(0) = r(1)*0.2387324146378430_dp
         omega(1) = r(2)*0.2387324146378430_dp
         omega(2) = r(3)*0.2387324146378430_dp
      case default
         xx(0) = 1.0_dp
         yy(0) = 1.0_dp
         zz(0) = 1.0_dp
         do i = 1, l
            xx(i) = xx(i - 1)*r(1)
            yy(i) = yy(i - 1)*r(2)
            zz(i) = zz(i - 1)*r(3)
         end do
         n = 0
         do i = l, 0, -1
            do j = l - i, 0, -1
               k = l - i - j
               omega(n) = xx(i)*yy(j)*zz(k)
               n = n + 1
            end do
         end do
         loc = cint_c2s_bra_sph(buf, 1, omega, l)
         call cint_s2c_bra_sph(buf, 1, omega, l)
      end select
   end subroutine ecp_ang_nuc_in_cart

   !> Integral of x^i y^j z^k over the unit sphere, divided by 4 pi
   !>
   !> Zero unless every exponent is even, which is what makes the angular sums
   !> above collapse: most terms never contribute and the parity tests in the
   !> callers exist to avoid forming them at all.
   pure function ecp_int_unit_xyz(i, j, k) result(v)
      integer, intent(in) :: i, j, k
      real(dp) :: v

      if (mod(i, 2) /= 0 .or. mod(j, 2) /= 0 .or. mod(k, 2) /= 0) then
         v = 0.0_dp
      else
         v = ecp_factorial2(i - 1)*ecp_factorial2(j - 1)*ecp_factorial2(k - 1) &
             /ecp_factorial2(i + j + k + 1)
      end if
   end function ecp_int_unit_xyz

   !> Binomial expansion factors for shifting a Cartesian shell to a new centre
   !>
   !> `facs` holds three (l+1)x(l+1) blocks, x then y then z, with
   !> `fac(i, j) = binom(i, j) * r^(i-j)`.
   pure subroutine ecp_cache_3dfac(facs, l, r)
      real(dp), intent(out) :: facs(0:)
      integer, intent(in) :: l
      real(dp), intent(in) :: r(3)

      real(dp) :: xx(0:15), yy(0:15), zz(0:15), bfac
      integer :: l1, i, j, off, ox, oy, oz

      l1 = l + 1
      ox = 0
      oy = l1*l1
      oz = 2*l1*l1

      xx(0) = 1.0_dp
      yy(0) = 1.0_dp
      zz(0) = 1.0_dp
      do i = 1, l
         xx(i) = xx(i - 1)*r(1)
         yy(i) = yy(i - 1)*r(2)
         zz(i) = zz(i - 1)*r(3)
      end do

      facs(0:3*l1*l1 - 1) = 0.0_dp
      do i = 0, l
         do j = 0, i
            bfac = ecp_binom(i, j)
            off = i*l1 + j
            facs(ox + off) = bfac*xx(i - j)
            facs(oy + off) = bfac*yy(i - j)
            facs(oz + off) = bfac*zz(i - j)
         end do
      end do
   end subroutine ecp_cache_3dfac

   !> Angular factors for the type-2 (projected) ECP term
   !>
   !> `facs` comes out shaped (li+1, nfi, dlc, dlambda) in the C's flat
   !> ordering, where dlc = 2*lc+1 and dlambda = li+lc+1.
   !>
   !> The parity test is the whole reason this is affordable. Only terms with
   !> (lc + a + b + c + lambda) even survive the angular integral, so half the
   !> lambda values are skipped rather than computed and discarded -- and the
   !> odd entries are zeroed explicitly, because the caller sums over the full
   !> range and would otherwise read whatever the workspace held.
   subroutine ecp_type2_facs_ang(facs, li, lc, ri)
      real(dp), intent(out) :: facs(0:)
      integer, intent(in) :: li, lc
      real(dp), intent(in) :: ri(3)

      real(dp) :: unitr(3)
      real(dp) :: omega_nuc(0:CART_CUM - 1)
      real(dp) :: buf(0:(14 + 1)*(14 + 2)/2 - 1)
      real(dp), allocatable :: omega(:), fac3d(:)
      integer :: li1, dlc, dlambda, dlclmb, nfi
      integer :: i, j, k, m, n, lmb, mi, loc
      integer :: need_even, need_odd, pu, pv, pw, pr, ps, pt
      integer :: base, po, pf, ox, oy, oz
      real(dp) :: acc, fac

      call unit_direction(ri, unitr)

      li1 = li + 1
      dlc = lc*2 + 1
      dlambda = li + lc + 1
      dlclmb = dlambda*dlc
      nfi = cint_len_cart(li)

      do i = 0, li + lc
         call ecp_ang_nuc_in_cart(omega_nuc(OFFSET_CART(i):), i, unitr)
      end do
      do i = 0, OFFSET_CART(li + lc + 1) - 1
         omega_nuc(i) = omega_nuc(i)*4.0_dp*acos(-1.0_dp)
      end do

      allocate (omega(0:li1*li1*li1*dlambda*dlc - 1))
      allocate (fac3d(0:3*li1*li1 - 1))

      do i = 0, li
         do j = 0, li - i
            do k = 0, li - i - j
               base = (i*li1*li1 + j*li1 + k)*dlclmb
               need_even = mod(lc + i + j + k, 2)

               po = base + need_even*dlc
               do lmb = need_even, li + lc, 2
                  do m = 0, cint_len_cart(lc) - 1
                     pv = CART_POW_Y(m)
                     pw = CART_POW_Z(m)
                     pu = lc - pv - pw
                     acc = 0.0_dp
                     do n = 0, cint_len_cart(lmb) - 1
                        ps = CART_POW_Y(n)
                        pt = CART_POW_Z(n)
                        pr = lmb - ps - pt
                        acc = acc + omega_nuc(OFFSET_CART(lmb) + n) &
                              *ecp_int_unit_xyz(i + pu + pr, j + pv + ps, k + pw + pt)
                     end do
                     buf(m) = acc
                  end do
                  select case (lc)
                  case (0)
                     omega(po) = buf(0)*0.282094791773878143_dp
                  case (1)
                     omega(po) = buf(0)*0.488602511902919921_dp
                     omega(po + 1) = buf(1)*0.488602511902919921_dp
                     omega(po + 2) = buf(2)*0.488602511902919921_dp
                  case default
                     loc = cint_c2s_bra_sph(omega(po:), 1, buf, lc)
                  end select
                  po = po + dlc*2
               end do

               need_odd = ieor(need_even, 1)
               po = base + need_odd*dlc
               do lmb = need_odd, li + lc, 2
                  do m = 0, dlc - 1
                     omega(po + m) = 0.0_dp
                  end do
                  po = po + dlc*2
               end do
            end do
         end do
      end do

      ox = 0
      oy = li1*li1
      oz = 2*li1*li1
      call ecp_cache_3dfac(fac3d, li, ri)

      facs(0:li1*nfi*dlclmb - 1) = 0.0_dp
      do mi = 0, nfi - 1
         ps = CART_POW_Y(mi)
         pt = CART_POW_Z(mi)
         pr = li - ps - pt
         do i = 0, pr
            do j = 0, ps
               do k = 0, pt
                  need_even = mod(lc + i + j + k, 2)
                  fac = fac3d(ox + pr*li1 + i)*fac3d(oy + ps*li1 + j)*fac3d(oz + pt*li1 + k)
                  base = (i*li1*li1 + j*li1 + k)*dlclmb
                  pf = ((i + j + k)*nfi + mi)*dlclmb
                  do m = 0, dlc - 1
                     do n = need_even, dlambda - 1, 2
                        facs(pf + m*dlambda + n) = facs(pf + m*dlambda + n) &
                                                   + fac*omega(base + n*dlc + m)
                     end do
                  end do
               end do
            end do
         end do
      end do

      deallocate (omega, fac3d)
   end subroutine ecp_type2_facs_ang

   !> Angular contraction of the type-1 (local) radial integrals
   !>
   !> `rad_all` arrives indexed by (total degree, lambda) and comes out
   !> contracted against the angular integral into `rad_ang`, shaped
   !> (lmax+1)^3 over the three Cartesian powers.
   subroutine ecp_type1_rad_ang(rad_ang, lmax, r, rad_all)
      real(dp), intent(out) :: rad_ang(0:)
      integer, intent(in) :: lmax
      real(dp), intent(in) :: r(3)
      real(dp), intent(in) :: rad_all(0:)

      real(dp) :: unitr(3)
      real(dp) :: omega_nuc(0:CART_CUM - 1)
      integer :: d1, d2, d3, i, j, k, n, lmb, need_even, pr, ps, pt
      real(dp) :: acc

      call unit_direction(r, unitr)

      do i = 0, lmax
         call ecp_ang_nuc_in_cart(omega_nuc(OFFSET_CART(i):), i, unitr)
      end do

      d1 = lmax + 1
      d2 = d1*d1
      d3 = d2*d1
      rad_ang(0:d3 - 1) = 0.0_dp

      do i = 0, lmax
         do j = 0, lmax - i
            do k = 0, lmax - i - j
               need_even = mod(i + j + k, 2)
               do lmb = need_even, lmax, 2
                  acc = 0.0_dp
                  do n = 0, cint_len_cart(lmb) - 1
                     ps = CART_POW_Y(n)
                     pt = CART_POW_Z(n)
                     pr = lmb - ps - pt
                     acc = acc + omega_nuc(OFFSET_CART(lmb) + n) &
                           *ecp_int_unit_xyz(i + pr, j + ps, k + pt)
                  end do
                  rad_ang(i*d2 + j*d1 + k) = rad_ang(i*d2 + j*d1 + k) &
                                             + rad_all((i + j + k)*d1 + lmb)*acc
               end do
            end do
         end do
      end do
   end subroutine ecp_type1_rad_ang

end module cint_ecp_ang
