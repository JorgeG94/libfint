! -*- fortran -*-
!
! Rys roots and weights by the Wheeler recursion.
!
! Included once per precision by cint_wheeler.F90.  The recursion turns a
! sequence of moments into the three-term coefficients of the orthogonal
! polynomials; the roots are then the eigenvalues of the resulting Jacobi
! matrix, and the weights come from the first components of its eigenvectors.
!
! Two moment generators, matching the C:
!
!   laguerre_moments  builds alpha/beta itself and works for any `lower`.
!   jacobi_moments    uses the fixed alpha/beta tables, and picks between
!                     Flocke's backward (Miller) recurrence for lower = 0 and
!                     the naive expansion in Boys functions otherwise.
!
! The eigenproblem is solved at double precision in both instantiations, as in
! the C: its extended paths also narrow to double before diagonalising.

   subroutine laguerre_moments(n, t, lower, alpha, beta, moments)
      integer,  intent(in)  :: n
      RKTYPE, intent(in)  :: t, lower
      RKTYPE, intent(out) :: alpha(0:), beta(0:), moments(0:)

      integer  :: i
      RKTYPE :: tt, t_inv, t2_inv, e0, l00, l01, fac0, fac1, l02
      RKTYPE :: lower2, l10, l11, l12, et

      tt = sqrt(t)
      t_inv = 0.5_rk / t
      t2_inv = 0.5_rk / (t * t)
      e0 = exp(-t) * t_inv
      l00 = 0.0_rk
      l01 = 1.0_rk

      alpha(0) = t_inv
      beta(0) = 0.0_rk

      if (lower == 0.0_rk) then
         moments(0) = sqrtpie4 / tt * erf(tt)
         moments(1) = -l01 * e0
         do i = 1, n - 2
            alpha(i) = real(i*4 + 1, rk) * t_inv
            beta(i)  = real(i * (i*2 - 1), rk) * t2_inv
            fac0 = real(i*4 - 1, rk) * t_inv
            fac1 = real((i - 1) * (i*2 - 1), rk) * t2_inv
            l02 = (1.0_rk - fac0) * l01 - fac1 * l00
            l00 = l01
            l01 = l02
            moments(i+1) = -l01 * e0
         end do
      else
         lower2 = lower * lower
         l10 = 0.0_rk
         l11 = 1.0_rk
         et = exp(-t * lower2) * lower * t_inv
         moments(0) = sqrtpie4 / tt * (erfc(lower * tt) - erfc(tt))
         moments(1) = l11 * et - l01 * e0
         do i = 1, n - 2
            alpha(i) = real(i*4 + 1, rk) * t_inv
            beta(i)  = real(i * (i*2 - 1), rk) * t2_inv
            fac0 = real(i*4 - 1, rk) * t_inv
            fac1 = real((i - 1) * (i*2 - 1), rk) * t2_inv
            l12 = (lower2 - fac0) * l11 - fac1 * l10
            l10 = l11
            l11 = l12
            l02 = (1.0_rk - fac0) * l01 - fac1 * l00
            l00 = l01
            l01 = l02
            moments(i+1) = l11 * et - l01 * e0
         end do
      end if
   end subroutine laguerre_moments

   ! Jacobi moments as an explicit combination of Boys functions.  The
   ! summation order comes from JACOBI_COEF_ORDER, which sorts the terms so
   ! that cancellation is as mild as it can be made -- so the order matters and
   ! is preserved.
   subroutine naive_jacobi_moments(n, t, lower, mus)
      integer,  intent(in)  :: n
      RKTYPE, intent(in)  :: t, lower
      RKTYPE, intent(out) :: mus(0:)

      integer  :: i, j, k, base
      RKTYPE :: s
      RKTYPE :: fmt(0:2*MXRYSROOTS-1)

      call fmt_erfc_like(fmt, t, lower, n - 1)

      do i = 0, n - 1
         base = i * (i + 1) / 2
         s = 0.0_rk
         do j = 0, i
            k = jacobi_coef_order(base + j)
            s = s + real(jacobi_coef(base + k), rk) * fmt(k)
         end do
         mus(i) = s
      end do
   end subroutine naive_jacobi_moments

   ! Flocke's recipe, JCP 131, 064107: run the recurrence backwards from an
   ! arbitrary seed well above the order wanted (Miller's algorithm), then fix
   ! the scale from the known zeroth moment.  Forward recurrence is unstable
   ! here; backward recurrence damps the unwanted solution instead.
   subroutine flocke_jacobi_moments(n, t, mus)
      integer,  intent(in)  :: n
      RKTYPE, intent(in)  :: t
      RKTYPE, intent(out) :: mus(0:)

      RKTYPE :: t_inv, mu0, mu1, mu2, rn, tt, norm
      integer  :: i

      if (t < smallx_limit) then
         call naive_jacobi_moments(n, t, 0.0_rk, mus)
         return
      end if

      t_inv = 0.5_rk / t
      mu1 = 1.0_rk    ! any non-zero seed will do
      mu2 = 0.0_rk
      mu0 = mu1

      do i = n - 1 + flocke_extra_order, n, -1
         rn = real(2*i + 3, rk) * t_inv + real(jacobi_rn_part2(i), rk)
         mu0 = (mu2 - rn * mu1) / real(jacobi_sn(i), rk)
         mu2 = mu1
         mu1 = mu0
      end do
      do i = n - 1, 0, -1
         rn = real(2*i + 3, rk) * t_inv + real(jacobi_rn_part2(i), rk)
         mu0 = (mu2 - rn * mu1) / real(jacobi_sn(i), rk)
         mus(i) = mu0
         mu2 = mu1
         mu1 = mu0
      end do

      tt = sqrt(t)
      norm = sqrtpie4 * erf(tt) / tt / mu0     ! fmt(0) / mu0
      do i = 0, n - 1
         mus(i) = mus(i) * norm
      end do
   end subroutine flocke_jacobi_moments

   ! Moments -> three-term recurrence coefficients.
   !
   ! The C carries a warning here worth keeping: when b of the three-term
   ! recurrence is small -- which happens as the lower limit of a short-range
   ! integral approaches 1, e.g. n=12, x=0.4, lower=0.93 -- sigma_kk decays
   ! like b^n and this recursion loses most of its digits.  The caller checks
   ! for it afterwards rather than the recursion avoiding it.
   subroutine wheeler_recursion(n, alpha, beta, moments, a, b)
      integer,  intent(in)  :: n
      RKTYPE, intent(in)  :: alpha(0:), beta(0:), moments(0:)
      RKTYPE, intent(out) :: a(0:), b(0:)

      integer  :: i, j, nc
      RKTYPE :: a0, b0, a1, b1
      ! Three rolling rows of the sigma table.  The C aliases pointers into one
      ! buffer and rotates them; separate arrays with an index rotation say the
      ! same thing without the aliasing.
      RKTYPE :: s0(0:2*MXRYSROOTS-1), sm(0:2*MXRYSROOTS-1), sk(0:2*MXRYSROOTS-1)
      RKTYPE :: tmp(0:2*MXRYSROOTS-1)

      a0 = alpha(0) + moments(1) / moments(0)
      b0 = 0.0_rk
      a(0) = a0
      b(0) = b0

      s0(0:2*n-1) = moments(0:2*n-1)
      sm = 0.0_rk
      sk = 0.0_rk

      do i = 1, n - 1
         nc = 2 * (n - i)
         do j = 0, nc - 1
            sk(j) = s0(2+j) - (a0 - alpha(i+j)) * s0(1+j) &
                    - b0 * sm(2+j) + beta(i+j) * s0(j)
         end do
         a1 = alpha(i) - s0(1) / s0(0) + sk(1) / sk(0)
         b1 = sk(0) / s0(0)
         a(i) = a1
         b(i) = b1
         a0 = a1
         b0 = b1
         ! rotate: sm <- s0, s0 <- sk, sk <- old sm
         tmp = sm
         sm = s0
         s0 = sk
         sk = tmp
      end do
   end subroutine wheeler_recursion

   ! Assemble the Jacobi matrix from the recurrence coefficients and solve it.
   !
   ! A non-positive b means the recursion above has broken down; the C truncates
   ! the rule there and zeroes the higher roots, treating their contribution as
   ! negligible, and this does the same.
   function rys_wheeler_partial(n_in, alpha, beta, moments, roots, weights) result(err)
      integer,  intent(in)    :: n_in
      RKTYPE, intent(in)    :: alpha(0:), beta(0:), moments(0:)
      real(dp), intent(out)   :: roots(0:), weights(0:)
      integer :: err

      RKTYPE :: a(0:MXRYSROOTS-1), b(0:MXRYSROOTS-1)
      real(dp) :: da(0:MXRYSROOTS-1), db(0:MXRYSROOTS-1)
      real(dp) :: c0(0:MXRYSROOTS*MXRYSROOTS-1)
      real(dp) :: mu0
      integer  :: i, k, n

      n = n_in
      mu0 = real(moments(0), dp)

      call wheeler_recursion(n, alpha, beta, moments, a, b)

      do i = 1, n - 1
         if (b(i) < 1.0e-14_rk) then
            if (b(i) < 0.0_rk) then
               ! Approximate the quadrature: higher-order roots go to zero.
               do k = i, n - 1
                  roots(k) = 0.0_dp
                  weights(k) = 0.0_dp
               end do
               n = i
               exit
            end if
         end if
         b(i) = sqrt(b(i))
      end do

      do i = 0, n - 1
         da(i) = real(a(i), dp)
      end do
      do i = 1, n - 1
         db(i-1) = real(b(i), dp)
      end do

      err = cint_diagonalize(n, da, db, roots, c0)

      do i = 0, n - 1
         roots(i) = roots(i) / (1.0_dp - roots(i))
         weights(i) = c0(i*n) * c0(i*n) * mu0
      end do
   end function rys_wheeler_partial

   function rys_laguerre(n, x, lower, roots, weights) result(err)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: x, lower
      real(dp), intent(out) :: roots(0:), weights(0:)
      integer :: err

      RKTYPE :: moments(0:2*MXRYSROOTS-1)
      RKTYPE :: alpha(0:2*MXRYSROOTS-1), beta(0:2*MXRYSROOTS-1)

      moments = 0.0_rk
      alpha = 0.0_rk
      beta = 0.0_rk
      call laguerre_moments(n*2, real(x, rk), real(lower, rk), alpha, beta, moments)
      err = rys_wheeler_partial(n, alpha, beta, moments, roots, weights)
   end function rys_laguerre

   function rys_jacobi(n, x, lower, roots, weights) result(err)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: x, lower
      real(dp), intent(out) :: roots(0:), weights(0:)
      integer :: err

      RKTYPE :: moments(0:2*MXRYSROOTS-1)

      moments = 0.0_rk
      if (lower == 0.0_dp) then
         call flocke_jacobi_moments(n*2, real(x, rk), moments)
      else
         call naive_jacobi_moments(n*2, real(x, rk), real(lower, rk), moments)
      end if
      err = rys_wheeler_partial(n, jacobi_alpha_rk, jacobi_beta_rk, moments, roots, weights)
   end function rys_jacobi
