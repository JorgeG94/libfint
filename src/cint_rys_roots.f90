!
! The Rys roots and weights entry points.
!
! This is the dispatch layer: two shortcuts at the ends of the x range, closed
! forms up to five roots, and above that a ladder of solvers chosen by root
! count and -- for the short-range kernel -- by the attenuation `lower`.
!
! The one structural difference from src/rys_roots.c is which solvers the arms
! point at.  The C has three precisions and so three sets of solver names; this
! has two, because the 80-bit and binary128 ladders collapse together
! (PORT_TO_FORTRAN.md 3.1).  Every arm the C sends to CINTlrys_* or CINTqrys_*
! goes to the quadruple module here, which is at least as accurate as either.
!
! Errors are returned, never fatal.  The C calls exit() unless KEEP_GOING is
! defined, which is not a reasonable thing for a library to do to its caller.
!
module cint_rys_roots
   use cint_const,        only: dp
   use cint_tab_roots_x0, only: POLY_SMALLX_R0, POLY_SMALLX_R1, &
                                POLY_SMALLX_W0, POLY_SMALLX_W1, &
                                POLY_LARGEX_RT, POLY_LARGEX_WW
   use cint_rys_root_n,   only: rys_root1, rys_root2, rys_root3, rys_root4, rys_root5
   use cint_schmidt_dp,   only: schmidt_dp => rys_schmidt
   use cint_schmidt_qp,   only: schmidt_qp => rys_schmidt
   use cint_wheeler_dp,   only: jacobi_dp => rys_jacobi, laguerre_dp => rys_laguerre
   use cint_wheeler_qp,   only: jacobi_qp => rys_jacobi, laguerre_qp => rys_laguerre
   implicit none
   private

   public :: cint_rys_roots_lr, cint_rys_roots_sr

   real(dp), parameter :: PIE4 = 0.78539816339744827900_dp
   real(dp), parameter :: SMALLX_LIMIT = 3.0e-7_dp

   ! Which solver an arm of the ladder uses.  Named rather than passed as
   ! procedure pointers: the set is closed and small, and a select case reads
   ! better than four abstract interfaces.
   integer, parameter :: S_SCHMIDT_DP  = 1
   integer, parameter :: S_SCHMIDT_QP  = 2
   integer, parameter :: S_JACOBI_DP   = 3
   integer, parameter :: S_JACOBI_QP   = 4
   integer, parameter :: S_LAGUERRE_DP = 5
   integer, parameter :: S_LAGUERRE_QP = 6

contains

   function solve(which, n, x, lower, u, w) result(err)
      integer,  intent(in)  :: which, n
      real(dp), intent(in)  :: x, lower
      real(dp), intent(out) :: u(0:), w(0:)
      integer :: err
      select case (which)
      case (S_SCHMIDT_DP);  err = schmidt_dp(n, x, lower, u, w)
      case (S_SCHMIDT_QP);  err = schmidt_qp(n, x, lower, u, w)
      case (S_JACOBI_DP);   err = jacobi_dp(n, x, lower, u, w)
      case (S_JACOBI_QP);   err = jacobi_qp(n, x, lower, u, w)
      case (S_LAGUERRE_DP); err = laguerre_dp(n, x, lower, u, w)
      case (S_LAGUERRE_QP); err = laguerre_qp(n, x, lower, u, w)
      case default;         err = 1
      end select
   end function solve

   ! Pick by x against a breakpoint, and fall back to the quadruple Schmidt if
   ! the chosen solver fails.
   function segment_solve(n, x, lower, u, w, breakpoint, f1, f2) result(err)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: x, lower, breakpoint
      integer,  intent(in)  :: f1, f2
      real(dp), intent(out) :: u(0:), w(0:)
      integer :: err

      if (x <= breakpoint) then
         err = solve(f1, n, x, lower, u, w)
      else
         err = solve(f2, n, x, lower, u, w)
      end if
      if (err /= 0) err = solve(S_SCHMIDT_QP, n, x, lower, u, w)
   end function segment_solve

   ! As above, but the arm is chosen by `lower` first and only then by x.
   ! Note the C does NOT retry here, unlike segment_solve.
   function segment_solve1(n, x, lower, u, w, lower_bp1, lower_bp2, breakpoint, &
                           f1, f2, f3) result(err)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: x, lower, lower_bp1, lower_bp2, breakpoint
      integer,  intent(in)  :: f1, f2, f3
      real(dp), intent(out) :: u(0:), w(0:)
      integer :: err

      if (lower < lower_bp1) then
         if (x <= breakpoint) then
            err = solve(f1, n, x, lower, u, w)
         else
            err = solve(f2, n, x, lower, u, w)
         end if
      else
         err = solve(f3, n, x, lower, u, w)
      end if
      ! lower_bp2 is accepted and unused, exactly as in the C, where it is
      ! passed by every caller and read by none.
      if (.false. .and. lower_bp2 > 0.0_dp) continue
   end function segment_solve1

   ! Long-range (unattenuated) Rys roots.
   function cint_rys_roots_lr(nroots, x, u, w) result(err)
      integer,  intent(in)  :: nroots
      real(dp), intent(in)  :: x
      real(dp), intent(out) :: u(0:), w(0:)
      integer :: err

      integer  :: off, i
      real(dp) :: rt, t

      err = 0
      if (x <= SMALLX_LIMIT) then
         off = nroots * (nroots - 1) / 2
         do i = 0, nroots - 1
            u(i) = POLY_SMALLX_R0(off+i) + POLY_SMALLX_R1(off+i) * x
            w(i) = POLY_SMALLX_W0(off+i) + POLY_SMALLX_W1(off+i) * x
         end do
         return
      else if (x >= 35.0_dp + nroots * 5.0_dp) then
         off = nroots * (nroots - 1) / 2
         t = sqrt(PIE4 / x)
         do i = 0, nroots - 1
            rt = POLY_LARGEX_RT(off+i)
            u(i) = rt / (x - rt)
            w(i) = POLY_LARGEX_WW(off+i) * t
         end do
         return
      end if

      select case (nroots)
      case (1);  err = rys_root1(x, u, w)
      case (2);  err = rys_root2(x, u, w)
      case (3);  err = rys_root3(x, u, w)
      case (4);  err = rys_root4(x, u, w)
      case (5);  err = rys_root5(x, u, w)
      case (6:7)
         err = segment_solve(nroots, x, 0.0_dp, u, w, 11.0_dp, S_JACOBI_DP, S_SCHMIDT_DP)
      case (8)
         err = segment_solve(nroots, x, 0.0_dp, u, w, 11.0_dp, S_JACOBI_DP, S_SCHMIDT_QP)
      case (9)
         err = segment_solve(nroots, x, 0.0_dp, u, w, 10.0_dp, S_JACOBI_QP, S_LAGUERRE_QP)
      case (10:11)
         err = segment_solve(nroots, x, 0.0_dp, u, w, 18.0_dp, S_JACOBI_QP, S_LAGUERRE_QP)
      case (12)
         err = segment_solve(nroots, x, 0.0_dp, u, w, 22.0_dp, S_JACOBI_QP, S_LAGUERRE_QP)
      case default
         err = segment_solve(nroots, x, 0.0_dp, u, w, 50.0_dp, S_JACOBI_QP, S_LAGUERRE_QP)
      end select
   end function cint_rys_roots_lr

   ! Short-range (erfc-attenuated) Rys roots.  `lower` is the attenuation
   ! ratio; the arms below are the C's, transcribed.
   function cint_rys_roots_sr(nroots, x, lower, u, w) result(err)
      integer,  intent(in)  :: nroots
      real(dp), intent(in)  :: x, lower
      real(dp), intent(out) :: u(0:), w(0:)
      integer :: err

      err = 0
      select case (nroots)
      case (1)
         err = schmidt_dp(nroots, x, lower, u, w)
      case (2)
         if (lower < 0.99_dp) then
            err = schmidt_dp(nroots, x, lower, u, w)
         else
            err = jacobi_qp(nroots, x, lower, u, w)
         end if
      case (3)
         if (lower < 0.93_dp) then
            err = schmidt_dp(nroots, x, lower, u, w)
         else if (lower < 0.97_dp) then
            err = segment_solve(nroots, x, lower, u, w, 10.0_dp, S_JACOBI_QP, S_LAGUERRE_QP)
         else
            err = jacobi_qp(nroots, x, lower, u, w)
         end if
      case (4)
         if (lower < 0.8_dp) then
            err = schmidt_dp(nroots, x, lower, u, w)
         else if (lower < 0.9_dp) then
            err = segment_solve(nroots, x, lower, u, w, 10.0_dp, S_JACOBI_QP, S_LAGUERRE_QP)
         else
            err = jacobi_qp(nroots, x, lower, u, w)
         end if
      case (5)
         if (lower < 0.4_dp) then
            err = segment_solve(nroots, x, lower, u, w, 50.0_dp, S_SCHMIDT_DP, S_LAGUERRE_QP)
         else if (lower < 0.8_dp) then
            err = segment_solve(nroots, x, lower, u, w, 10.0_dp, S_JACOBI_QP, S_LAGUERRE_QP)
         else
            err = jacobi_qp(nroots, x, lower, u, w)
         end if
      case (6)
         if (lower < 0.25_dp) then
            err = segment_solve(nroots, x, lower, u, w, 60.0_dp, S_SCHMIDT_DP, S_LAGUERRE_QP)
         else if (lower < 0.8_dp) then
            err = segment_solve(nroots, x, lower, u, w, 10.0_dp, S_JACOBI_QP, S_LAGUERRE_QP)
         else
            err = jacobi_qp(nroots, x, lower, u, w)
         end if
      case (7)
         err = segment_solve1(nroots, x, lower, u, w, 0.5_dp, 1.0_dp, 60.0_dp, &
                              S_JACOBI_QP, S_LAGUERRE_QP, S_JACOBI_QP)
      case (8:12)
         err = segment_solve1(nroots, x, lower, u, w, 0.15_dp, 1.0_dp, 60.0_dp, &
                              S_JACOBI_QP, S_LAGUERRE_QP, S_JACOBI_QP)
      case (13:14)
         err = segment_solve1(nroots, x, lower, u, w, 0.25_dp, 1.0_dp, 60.0_dp, &
                              S_JACOBI_QP, S_LAGUERRE_QP, S_JACOBI_QP)
      case (15:16)
         err = segment_solve1(nroots, x, lower, u, w, 0.25_dp, 0.75_dp, 60.0_dp, &
                              S_JACOBI_QP, S_LAGUERRE_QP, S_JACOBI_QP)
      case (17)
         err = segment_solve1(nroots, x, lower, u, w, 0.25_dp, 0.65_dp, 60.0_dp, &
                              S_JACOBI_QP, S_LAGUERRE_QP, S_JACOBI_QP)
      case (18)
         err = segment_solve1(nroots, x, lower, u, w, 0.15_dp, 0.65_dp, 60.0_dp, &
                              S_JACOBI_QP, S_LAGUERRE_QP, S_JACOBI_QP)
      case (19)
         err = segment_solve1(nroots, x, lower, u, w, 0.15_dp, 0.55_dp, 60.0_dp, &
                              S_JACOBI_QP, S_LAGUERRE_QP, S_JACOBI_QP)
      case (20:21)
         err = segment_solve1(nroots, x, lower, u, w, 0.25_dp, 0.45_dp, 60.0_dp, &
                              S_JACOBI_QP, S_LAGUERRE_QP, S_JACOBI_QP)
      case (22:24)
         err = segment_solve1(nroots, x, lower, u, w, 0.25_dp, 0.35_dp, 60.0_dp, &
                              S_JACOBI_QP, S_LAGUERRE_QP, S_JACOBI_QP)
      case default
         err = 1
      end select
   end function cint_rys_roots_sr

end module cint_rys_roots
