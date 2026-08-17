! -*- fortran -*-
!
! The Boys function F_m(t) and its short-range (erfc-attenuated) form.
!
! This file is the ALGORITHM ONLY.  It is included once per precision by
! cint_fmt.F90, into a module that has already bound `rk` to a kind and `sml`
! to that kind's convergence tolerance.  Fortran has no generics, and the C
! solves the same problem by writing the routines out three times -- once for
! double, once for long double, once for __float128.  Keeping one copy and
! including it twice is the difference between one algorithm and three
! transcriptions that can drift apart.
!
! Two rather than three because the long double and __float128 ladders both
! become binary128 here; see PORT_TO_FORTRAN.md 3.1.
!
! Indices are 0-based to match the C: f(0:m) holds F_0 .. F_m.

   ! F_m(t) by the ascending series.  Used below the turnover point, where
   ! the recursion downward from F_m is stable and the closed form is not.
   subroutine fmt1_gamma_inc_like(f, t, m)
      real(rk), intent(out) :: f(0:)
      real(rk), intent(in)  :: t
      integer,  intent(in)  :: m

      integer  :: i
      real(rk) :: b, bi, e, x, s, tol

      b = real(m, rk) + 0.5_rk
      e = 0.5_rk * exp(-t)
      x = e
      s = e
      tol = sml * e
      bi = b + 1.0_rk
      do while (x > tol)
         x = x * (t / bi)
         s = s + x
         bi = bi + 1.0_rk
      end do
      f(m) = s / b
      do i = m, 1, -1
         b = b - 1.0_rk
         f(i-1) = (e + t * f(i)) / b
      end do
   end subroutine fmt1_gamma_inc_like

   subroutine gamma_inc_like(f, t, m)
      real(rk), intent(out) :: f(0:)
      real(rk), intent(in)  :: t
      integer,  intent(in)  :: m

      integer  :: i
      real(rk) :: tt, e, b

      if (t == 0.0_rk) then
         f(0) = 1.0_rk
         do i = 1, m
            f(i) = 1.0_rk / real(2*i + 1, rk)
         end do
      else if (t < turnover_point(m)) then
         call fmt1_gamma_inc_like(f, t, m)
      else
         tt = sqrt(t)
         f(0) = sqrtpie4 / tt * erf(tt)
         e = exp(-t)
         b = 0.5_rk / t
         do i = 1, m
            f(i) = b * (real(2*i - 1, rk) * f(i-1) - e)
         end do
      end if
   end subroutine gamma_inc_like

   ! Series form of the attenuated integral.  The subtraction (x - x1) is the
   ! delicate part: both terms are positive and nearly equal as `lower`
   ! approaches 1, which is why the C documents this whole path as unstable
   ! for x*lower**2 above about 40.
   subroutine fmt1_erfc_like(f, t, lower, m)
      real(rk), intent(out) :: f(0:)
      real(rk), intent(in)  :: t, lower
      integer,  intent(in)  :: m

      integer  :: i
      real(rk) :: lower2, b, bi, e, e1, x, x1, s, div, delta, tol, val

      lower2 = lower * lower
      b = real(m, rk) + 0.5_rk
      e = 0.5_rk * exp(-t)
      e1 = 0.5_rk * exp(-t * lower2) * lower
      e1 = e1 * lower2 ** m
      x = e
      x1 = e1
      s = e - e1
      div = 1.0_rk
      delta = s
      tol = sml * abs(delta)
      bi = b + 1.0_rk
      do while (abs(delta) > tol)
         div = div * (t / bi)
         x1 = x1 * lower2
         delta = (x - x1) * div
         s = s + delta
         bi = bi + 1.0_rk
      end do
      val = s / b
      f(m) = val
      do i = m, 1, -1
         b = b - 1.0_rk
         e1 = e1 / lower2
         val = (e - e1 + t * val) / b
         f(i-1) = val
      end do
   end subroutine fmt1_erfc_like

   subroutine fmt_erfc_like(f, t, lower, m)
      real(rk), intent(out) :: f(0:)
      real(rk), intent(in)  :: t, lower
      integer,  intent(in)  :: m

      integer  :: i
      real(rk) :: lower2, e1, tt, val, e, b

      if (lower == 0.0_rk) then
         call gamma_inc_like(f, t, m)
         return
      end if

      lower2 = lower * lower
      ! F(m) < .5*sqrt(pi/t) * erfc(low*tt); below this the answer is zero to
      ! working precision and the series would only manufacture noise.
      if (t * lower2 > erfc_bound) then
         do i = 0, m
            f(i) = 0.0_rk
         end do
         return
      end if

      if (t == 0.0_rk) then
         e1 = lower
         f(0) = 1.0_rk - e1
         do i = 1, m
            e1 = e1 * lower2
            f(i) = (1.0_rk - e1) / real(2*i + 1, rk)
         end do
      else if (t < turnover_point(m)) then
         call fmt1_erfc_like(f, t, lower, m)
      else
         tt = sqrt(t)
         ! erfc(a) - erfc(b) is more accurate here than erf(b) - erf(a)
         val = sqrtpie4 / tt * (erfc(lower * tt) - erfc(tt))
         f(0) = val
         e = exp(-t)
         e1 = exp(-t * lower2) * lower
         b = 0.5_rk / t
         do i = 0, m - 1
            val = b * (real(2*i + 1, rk) * val - e + e1)
            e1 = e1 * lower2
            f(i+1) = val
         end do
      end if
   end subroutine fmt_erfc_like
