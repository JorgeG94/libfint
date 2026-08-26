!
! Radial parts of the ECP integrals.
!
! Translated from PySCF's pyscf/lib/gto/nr_ecp.c.  See cint_ecp_num for why
! the upstream here is PySCF rather than libcint.
!
! **This is where bit-identity with PySCF stops, deliberately.**  The C calls
! `ECPsph_ine_opt` from both radial routines below, which is a seventh-order
! Taylor interpolation over a 400-entry table -- about 3,600 lines of
! coefficients whose only purpose is to be faster than the series.  This port
! calls the series instead, and the difference is worth stating precisely
! rather than leaving as an asterisk:
!
!   * over z in [1e-7, 16) the two agree to 4.4e-16 absolute, worst case,
!     on values of order one.  That is two or three ULP, so anything built
!     on this matches PySCF to about 1e-15 relative and not to the last bit.
!   * the series is the more accurate of the two.  The interpolation loses
!     relative accuracy badly for tiny values -- at z = 1.6e-4 and l = 4 it
!     returns -8.8e-19 where the true value is +7.1e-19, the wrong sign -- and
!     an ECP integral only escapes that because such terms are negligible
!     next to the ones that matter.
!   * at exactly z = 16 the C reads past its table.  `ECPsph_ine_opt`
!     delegates to the series for `z > 16`, so z = 16 itself takes the table
!     branch and indexes entry floor(16/0.04) = 400 of a 400-entry table.
!     Reproducing that is not something a translation should do.
!
! So the tables stay out and the bar below this point is a tolerance rather
! than an equality.  cint_ecp_num, which is upstream of this choice, is still
! checked bit for bit.
!
module cint_ecp_rad
   use cint_const, only: dp
   use cint_bas, only: ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, BAS_SLOTS
   use cint_ecp_num, only: ecp_sph_ine, ECP_SIM_ZERO, ECP_EXPCUTOFF, ECP_CUTOFF, &
                           RADI_POWER
   use cint_ecp_ang, only: ecp_cache_3dfac, CART_POW_Y, CART_POW_Z
   use cint_bas, only: cint_len_cart
   implicit none
   private

   public :: ecp_rad_part
   public :: ecp_type2_facs_rad
   public :: ecp_type1_static_facs
   public :: ecp_type1_rad_part

contains

   !> The ECP's radial potential U(r) summed over one group of ECP shells
   !>
   !> Returns the number of grid points that were still non-negligible, which
   !> the callers use to shorten every loop that follows.  That truncation is
   !> the main saving in the whole evaluation: an ECP is short ranged, so on a
   !> 2047-point grid most of the tail contributes nothing.
   !>
   !> The break test needs two consecutive small values, not one. A single
   !> small value is not evidence of decay -- a contracted potential can pass
   !> through zero -- and stopping there would truncate a live integrand.
   subroutine ecp_rad_part(ur, nrs_max, rs, rs_off, nrs, inc, ecpsh0, ecpsh1, &
                           ecpbas, env)
      real(dp), intent(out) :: ur(0:)
      integer, intent(out) :: nrs_max
      real(dp), intent(in) :: rs(0:)
      integer, intent(in) :: rs_off, nrs, inc, ecpsh0, ecpsh1
      integer, intent(in) :: ecpbas(0:)
      real(dp), intent(in) :: env(0:)

      real(dp), allocatable :: ubuf(:), r2(:)
      real(dp) :: s
      integer :: ish, i, kp, n, npk, nrs_now, pexp, pcoeff, rpow

      allocate (ubuf(0:nrs - 1), r2(0:nrs - 1))
      nrs_max = 0

      do i = 0, nrs - 1
         r2(i) = rs(rs_off + i*inc)*rs(rs_off + i*inc)
         ur(i) = 0.0_dp
      end do

      do ish = ecpsh0, ecpsh1 - 1
         npk = ecpbas(ish*BAS_SLOTS + NPRIM_OF)
         pexp = ecpbas(ish*BAS_SLOTS + PTR_EXP)
         pcoeff = ecpbas(ish*BAS_SLOTS + PTR_COEFF)

         nrs_now = nrs
         do i = 0, nrs - 1
            s = env(pcoeff)*exp(-env(pexp)*r2(i))
            do kp = 1, npk - 1
               s = s + env(pcoeff + kp)*exp(-env(pexp + kp)*r2(i))
            end do
            ubuf(i) = s
            if (i > 2 .and. abs(ubuf(i)) < ECP_SIM_ZERO &
                .and. abs(ubuf(i - 1)) < ECP_SIM_ZERO) then
               nrs_now = i
               exit
            end if
         end do

         nrs_max = max(nrs_max, nrs_now)

         rpow = ecpbas(ish*BAS_SLOTS + RADI_POWER)
         select case (rpow)
         case (1)
            do i = 0, nrs_now - 1
               ubuf(i) = ubuf(i)*rs(rs_off + i*inc)
            end do
         case (2)
            do i = 0, nrs_now - 1
               ubuf(i) = ubuf(i)*r2(i)
            end do
         case (3)
            do i = 0, nrs_now - 1
               ubuf(i) = ubuf(i)*(r2(i)*rs(rs_off + i*inc))
            end do
         case default
            do i = 0, nrs_now - 1
               do n = 1, rpow
                  ubuf(i) = ubuf(i)*rs(rs_off + i*inc)
               end do
            end do
         end select

         do i = 0, nrs_now - 1
            ur(i) = ur(i) + ubuf(i)
         end do
      end do

      deallocate (ubuf, r2)
   end subroutine ecp_rad_part

   !> Radial factors for the type-2 (projected) term, contracted over primitives
   !>
   !> `facs` comes out shaped (nrs, li+lc+1, nc): for each contracted
   !> function, each grid point carries the whole lambda ladder.
   subroutine ecp_type2_facs_rad(facs, ish, lc, rca, rs, nrs, inc, bas, env)
      real(dp), intent(out) :: facs(0:)
      integer, intent(in) :: ish, lc, nrs, inc
      real(dp), intent(in) :: rca
      real(dp), intent(in) :: rs(0:)
      integer, intent(in) :: bas(0:)
      real(dp), intent(in) :: env(0:)

      integer :: li, np, nc, lilc1, ip, i, j, ic, m, pexp, pcoeff
      real(dp) :: ka, t1, ar2
      real(dp), allocatable :: r2(:), buf(:), bess(:)
      real(dp) :: acc

      if (nrs == 0) return

      li = bas(ish*BAS_SLOTS + ANG_OF)
      np = bas(ish*BAS_SLOTS + NPRIM_OF)
      nc = bas(ish*BAS_SLOTS + NCTR_OF)
      pexp = bas(ish*BAS_SLOTS + PTR_EXP)
      pcoeff = bas(ish*BAS_SLOTS + PTR_COEFF)
      lilc1 = li + lc + 1

      allocate (r2(0:nrs - 1), buf(0:np*nrs*lilc1 - 1), bess(0:lilc1 - 1))

      do i = 0, nrs - 1
         t1 = rs(i*inc) - rca
         r2(i) = t1*t1
      end do

      do ip = 0, np - 1
         ka = 2.0_dp*env(pexp + ip)*rca
         do i = 0, nrs - 1
            ar2 = env(pexp + ip)*r2(i)
            ! The + 6 is the C's, and is about the largest value the Bessel
            ! routine can still return usefully rather than about exp itself.
            if (ar2 > ECP_EXPCUTOFF + 6.0_dp) then
               buf((ip*nrs + i)*lilc1:(ip*nrs + i)*lilc1 + li + lc) = 0.0_dp
            else
               t1 = exp(-ar2)
               call ecp_sph_ine(bess, li + lc, ka*rs(i*inc))
               do j = 0, li + lc
                  buf((ip*nrs + i)*lilc1 + j) = bess(j)*t1
               end do
            end if
         end do
      end do

      ! The C's dgemm over (nrs*lilc1) x nc x np.  Written out because the
      ! left operand is the buffer just built and the right is nc columns of
      ! contraction coefficients; at these sizes a call would cost more than
      ! the loop.
      m = nrs*lilc1
      do ic = 0, nc - 1
         do i = 0, m - 1
            acc = 0.0_dp
            do ip = 0, np - 1
               acc = acc + buf(ip*m + i)*env(pcoeff + ic*np + ip)
            end do
            facs(ic*m + i) = acc
         end do
      end do

      deallocate (r2, buf, bess)
   end subroutine ecp_type2_facs_rad

   !> Centre-shift factors for the type-1 term, one block per Cartesian component
   subroutine ecp_type1_static_facs(facs, li, ri)
      real(dp), intent(out) :: facs(0:)
      integer, intent(in) :: li
      real(dp), intent(in) :: ri(3)

      real(dp), allocatable :: fac3d(:)
      integer :: d1, d2, d3, mi, i, j, k, px, py, pz, ox, oy, oz

      d1 = li + 1
      d2 = d1*d1
      d3 = d2*d1
      ox = 0
      oy = d1*d1
      oz = 2*d1*d1

      allocate (fac3d(0:3*d1*d1 - 1))
      call ecp_cache_3dfac(fac3d, li, ri)

      facs(0:cint_len_cart(li)*d3 - 1) = 0.0_dp
      do mi = 0, cint_len_cart(li) - 1
         py = CART_POW_Y(mi)
         pz = CART_POW_Z(mi)
         px = li - py - pz
         do i = 0, px
            do j = 0, py
               do k = 0, pz
                  facs(mi*d3 + i*d2 + j*d1 + k) = &
                     fac3d(ox + px*d1 + i)*fac3d(oy + py*d1 + j)*fac3d(oz + pz*d1 + k)
               end do
            end do
         end do
      end do

      deallocate (fac3d)
   end subroutine ecp_type1_static_facs

   !> Radial integrals for the type-1 (local) term
   !>
   !> Accumulates into `rad_all`, which the caller has zeroed: this is called
   !> once per primitive pair and the sum over them is what builds the result.
   !>
   !> The parity stride in the innermost loop is not an optimisation. Terms
   !> where (lab + i) is odd integrate to zero over the sphere, so stepping by
   !> two skips exactly the ones the angular factor would multiply by nothing.
   subroutine ecp_type1_rad_part(rad_all, lmax, k, aij, ur, rs, nrs, inc)
      real(dp), intent(inout) :: rad_all(0:)
      integer, intent(in) :: lmax, nrs, inc
      real(dp), intent(in) :: k, aij
      real(dp), intent(in) :: ur(0:), rs(0:)

      integer :: lmax1, lab, i, n
      real(dp) :: kaij, fac, tmp, s
      real(dp), allocatable :: rur(:), bval(:), bess(:)

      if (nrs == 0) return

      lmax1 = lmax + 1
      allocate (rur(0:nrs - 1), bval(0:nrs*lmax1 - 1), bess(0:lmax))

      kaij = k/(2.0_dp*aij)
      fac = kaij*kaij*aij
      do n = 0, nrs - 1
         tmp = rs(n*inc) - kaij
         tmp = fac - aij*tmp*tmp
         ! Three ways to be negligible, and the first two are about the
         ! arithmetic rather than the physics: exp(tmp) would overflow above
         ! CUTOFF, and underflow below. The C's comment is that such points
         ! are remote functions and dropping them costs no accuracy.
         if (ur(n) == 0.0_dp .or. tmp > ECP_CUTOFF &
             .or. tmp < -(ECP_EXPCUTOFF + 6.0_dp + 30.0_dp)) then
            rur(n) = 0.0_dp
            bval(n*lmax1:n*lmax1 + lmax) = 0.0_dp
         else
            rur(n) = ur(n)*exp(tmp)
            call ecp_sph_ine(bess, lmax, k*rs(n*inc))
            bval(n*lmax1:n*lmax1 + lmax) = bess(0:lmax)
         end if
      end do

      do lab = 0, lmax
         if (lab > 0) then
            do n = 0, nrs - 1
               rur(n) = rur(n)*rs(n*inc)
            end do
         end if
         do i = mod(lab, 2), lmax, 2
            s = rad_all(lab*lmax1 + i)
            do n = 0, nrs - 1
               s = s + rur(n)*bval(n*lmax1 + i)
            end do
            rad_all(lab*lmax1 + i) = s
         end do
      end do

      deallocate (rur, bval, bess)
   end subroutine ecp_type1_rad_part

end module cint_ecp_rad
