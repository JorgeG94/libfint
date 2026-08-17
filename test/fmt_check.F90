!
! Validation for the Fortran Boys function (deliverable D3, fmt layer).
!
! Two independent checks, because neither alone is sufficient.
!
! 1. TRANSCRIPTION.  The double-precision path is compared against a dump from
!    the C's own gamma_inc_like / fmt_erfc_like.  This asks only "did the port
!    copy the algorithm faithfully", and the answer should be that every value
!    agrees to the last bit.  It says nothing about whether the algorithm is
!    any good.
!
! 2. CORRECTNESS.  The quadruple path is checked against mathematics rather
!    than against the C, for the good reason that the C is not accurate enough
!    to serve as a reference here.  At t = 1e-12, lower = 0.75, m = 1 the true
!    value of F_1 is (1 - lower^3)/3 = 0.19270833...; the C's double returns
!    -11.12 and its long double returns 0.2142, while this code returns
!    0.1927083333332.  Comparing against either C path would flag correct
!    answers as failures.
!
!    So instead: the exact closed forms at t = 0, and the recurrence
!
!        (2m-1) F_{m-1} = e^-t - lower^(2m-1) e^(-t lower^2) + 2t F_m
!
!    which every F_m must satisfy regardless of how it was computed.  A wrong
!    implementation has to be wrong in a very particular, self-consistent way
!    to survive that.
!
program fmt_check
   use cint_const,  only: dp, qp
   use cint_fmt_dp, only: gamma_inc_like, fmt_erfc_like
   use cint_fmt_qp, only: q_gamma_inc_like => gamma_inc_like, &
                          q_fmt_erfc_like  => fmt_erfc_like
   implicit none

   integer,  parameter :: MMAX = 20
   real(dp), parameter :: TS(22) = [ &
      0.0_dp, 1e-12_dp, 1e-6_dp, 1e-3_dp, 0.01_dp, 0.1_dp, 0.5_dp, &
      0.86602540378_dp, 0.87_dp, 1.0_dp, 1.295010032056_dp, 1.3_dp, &
      2.0_dp, 3.28_dp, 5.0_dp, 8.0_dp, 12.0_dp, 20.0_dp, 40.0_dp, &
      100.0_dp, 200.0_dp, 500.0_dp ]
   real(dp), parameter :: LS(8) = [ &
      0.0_dp, 0.1_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.9_dp, 0.95_dp, 0.99_dp ]

   logical :: ok
   integer :: nfail

   nfail = 0

   call check_transcription(ok)
   if (.not. ok) nfail = nfail + 1

   call check_recurrence(ok)
   if (.not. ok) nfail = nfail + 1

   call check_exact_at_zero(ok)
   if (.not. ok) nfail = nfail + 1

   print '(A)', ""
   if (nfail == 0) then
      print '(A)', "  RESULT: PASS"
   else
      print '(A,I0,A)', "  RESULT: FAIL (", nfail, " checks failed)"
      stop 1
   end if

contains

   ! ---- 1. does the double path reproduce the C exactly? ----------------
   subroutine check_transcription(ok)
      logical, intent(out) :: ok
      real(dp) :: f(0:MMAX+1), c(0:MMAX+1)
      integer  :: m, i, j, k, u, n, ios
      real(dp) :: rel, worst
      character(len=256) :: path

      ok = .true.
      call get_command_argument(1, path)
      open(newunit=u, file=trim(path), access='stream', form='unformatted', &
           status='old', iostat=ios)
      if (ios /= 0) then
         print '(A,A)', "  cannot open C reference dump: ", trim(path)
         ok = .false.
         return
      end if

      worst = 0.0_dp
      n = 0
      do m = 0, MMAX
         do i = 1, size(TS)
            read(u) c(0:m)
            call gamma_inc_like(f, TS(i), m)
            do k = 0, m
               rel = reldiff(f(k), c(k))
               worst = max(worst, rel)
               n = n + 1
            end do
            do j = 1, size(LS)
               read(u) c(0:m)
               call fmt_erfc_like(f, TS(i), LS(j), m)
               do k = 0, m
                  rel = reldiff(f(k), c(k))
                  worst = max(worst, rel)
                  n = n + 1
               end do
            end do
         end do
      end do
      close(u)

      print '(A)',            "  [1] transcription: dp against the C's own double"
      print '(A,I0)',         "      values compared : ", n
      print '(A,ES10.2)',     "      worst rel diff  : ", worst
      ! Bitwise is the expectation: same algorithm, same order of operations,
      ! same precision.  Anything above zero means the port drifted.
      if (worst /= 0.0_dp) then
         print '(A)',         "      EXPECTED 0.0 -- the port has drifted from the C"
         ok = .false.
      end if
   end subroutine check_transcription

   ! ---- 2. does the quad path satisfy the Boys recurrence? --------------
   subroutine check_recurrence(ok)
      logical, intent(out) :: ok
      real(qp) :: f(0:MMAX+1), t, l, e, e1, lhs, rhs, scale
      integer  :: m, i, j
      real(qp) :: rel, worst_g, worst_e
      integer  :: nbad

      ok = .true.
      worst_g = 0.0_qp
      worst_e = 0.0_qp
      nbad = 0

      do i = 1, size(TS)
         t = real(TS(i), qp)
         if (t == 0.0_qp) cycle
         ! unattenuated: (2m-1) F_{m-1} = e^-t + 2t F_m
         call q_gamma_inc_like(f, t, MMAX)
         e = exp(-t)
         do m = 1, MMAX
            lhs = real(2*m - 1, qp) * f(m-1)
            rhs = e + 2.0_qp * t * f(m)
            scale = max(abs(lhs), abs(rhs))
            if (scale > 0.0_qp) then
               rel = abs(lhs - rhs) / scale
               worst_g = max(worst_g, rel)
               if (rel > 1.0e-28_qp) nbad = nbad + 1
            end if
         end do

         ! attenuated: (2m-1) F_{m-1} = e^-t - l^(2m-1) e^(-t l^2) + 2t F_m
         do j = 1, size(LS)
            l = real(LS(j), qp)
            if (l == 0.0_qp) cycle
            call q_fmt_erfc_like(f, t, l, MMAX)
            if (all(f(0:MMAX) == 0.0_qp)) cycle   ! the erfc_bound cutoff
            e1 = exp(-t * l * l)
            do m = 1, MMAX
               lhs = real(2*m - 1, qp) * f(m-1)
               rhs = e - l**(2*m - 1) * e1 + 2.0_qp * t * f(m)
               scale = max(abs(lhs), abs(rhs))
               if (scale > 0.0_qp) then
                  rel = abs(lhs - rhs) / scale
                  worst_e = max(worst_e, rel)
                  if (rel > 1.0e-24_qp) nbad = nbad + 1
               end if
            end do
         end do
      end do

      print '(A)',        "  [2] correctness: qp against the Boys recurrence"
      print '(A,ES10.2)', "      worst residual, unattenuated : ", real(worst_g, dp)
      print '(A,ES10.2)', "      worst residual, attenuated   : ", real(worst_e, dp)
      if (nbad > 0) then
         print '(A,I0)',  "      residuals beyond tolerance   : ", nbad
         ok = .false.
      end if
   end subroutine check_recurrence

   ! ---- 3. exact closed forms at t = 0 ----------------------------------
   subroutine check_exact_at_zero(ok)
      logical, intent(out) :: ok
      real(qp) :: f(0:MMAX+1), l, want
      real(dp) :: rel, worst
      integer  :: m, j, nbad

      ok = .true.
      worst = 0.0_dp
      nbad = 0

      ! F_m(0) = 1/(2m+1)
      call q_gamma_inc_like(f, 0.0_qp, MMAX)
      do m = 0, MMAX
         want = 1.0_qp / real(2*m + 1, qp)
         rel = real(abs(f(m) - want) / want, dp)
         worst = max(worst, rel)
         if (rel > 1.0e-30_dp) nbad = nbad + 1
      end do

      ! attenuated at t = 0: F_m = (1 - l^(2m+1)) / (2m+1)
      do j = 1, size(LS)
         l = real(LS(j), qp)
         if (l == 0.0_qp) cycle
         call q_fmt_erfc_like(f, 0.0_qp, l, MMAX)
         do m = 0, MMAX
            want = (1.0_qp - l**(2*m + 1)) / real(2*m + 1, qp)
            rel = real(abs(f(m) - want) / want, dp)
            worst = max(worst, rel)
            if (rel > 1.0e-30_dp) nbad = nbad + 1
         end do
      end do

      print '(A)',        "  [3] correctness: qp against the exact values at t = 0"
      print '(A,ES10.2)', "      worst rel diff  : ", worst
      if (nbad > 0) then
         print '(A,I0)',  "      beyond tolerance: ", nbad
         ok = .false.
      end if
   end subroutine check_exact_at_zero

   pure real(dp) function reldiff(a, b)
      real(dp), intent(in) :: a, b
      if (b /= 0.0_dp) then
         reldiff = abs(a - b) / abs(b)
      else
         reldiff = abs(a - b)
      end if
   end function reldiff

end program fmt_check
