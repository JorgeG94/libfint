!
! Boys function for the Rys quadrature, at both working precisions.
!
! The algorithm lives once, in cint_fmt_body.inc, and is included twice --
! into a double module and a quadruple one.  Each include sees a different
! `rk` and a different convergence tolerance and produces its own copy of the
! four routines; because the copies live in separate modules, they can share
! the same names without any preprocessor token pasting.
!
! Precision map against the C (see PORT_TO_FORTRAN.md 3.1):
!
!   C gamma_inc_like    double        ->  cint_fmt_dp
!   C lgamma_inc_like   long double   ->  cint_fmt_qp   (80-bit -> binary128)
!   C qgamma_inc_like   __float128    ->  cint_fmt_qp
!
! The middle row is a promotion, so the Fortran is more accurate than the C
! wherever the long double ladder is used, and its tolerance is the C's
! SML_FLOAT128 rather than SML_FLOAT80 for the same reason.
!
module cint_fmt_tab
   use cint_const, only: dp
   implicit none
   public

   ! Below TURNOVER_POINT(m) the ascending series wins; above it the closed
   ! form with erf() is both faster and better conditioned.
   real(dp), parameter :: turnover_point(0:39) = [ &
      0.0_dp, 0.0_dp, 0.866025403784_dp, 1.295010032056_dp, 1.705493613097_dp, 2.106432965305_dp, &
      2.501471934009_dp, 2.892473348218_dp, 3.280525047072_dp, 3.666320693281_dp, 4.05033123037_dp, &
      4.432891808508_dp, 4.814249856864_dp, 5.194593501454_dp, 5.574069276051_dp, 5.952793645111_dp, &
      6.330860773135_dp, 6.708347923415_dp, 7.08531930745_dp, 7.461828891625_dp, 7.837922483937_dp, &
      8.213639312398_dp, 8.589013237349_dp, 8.964073695432_dp, 9.338846443746_dp, 9.713354153046_dp, &
      10.08761688545_dp, 10.46165248270_dp, 10.83547688448_dp, 11.20910439128_dp, 11.58254788331_dp, &
      11.95581900374_dp, 12.32892831326_dp, 12.70188542111_dp, 13.07469909673_dp, 13.44737736550_dp, &
      13.81992759110_dp, 14.19235654675_dp, 14.56467047710_dp, 14.93687515212_dp ]
end module cint_fmt_tab

module cint_fmt_dp
   ! `rk => dp` comes FIRST here, and the order is not cosmetic.  LFortran
   ! 0.64 keeps only the LAST entry of an only-list that names a given
   ! symbol, so `only: dp, rk => dp` leaves `dp` itself undeclared in this
   ! module -- and reports it against cint_quad.f90, which does not have the
   ! problem, so the message points nowhere near the cause.  Written this way
   ! round both names survive.  Every other compiler reads the two orders
   ! identically.  The same import appears five more times -- below, and
   ! twice each in cint_schmidt.F90 and cint_wheeler.F90 -- each with a
   ! one-line pointer back to this note.
   use cint_const,  only: rk => dp, dp
#undef RKTYPE
#undef TO_RK
#undef TO_DP
#define RKTYPE real(rk)
#define TO_RK(x) real(x, rk)
#define TO_DP(x) real(x, dp)
   use cint_fmt_tab, only: turnover_point
   implicit none
   private
   public :: gamma_inc_like, fmt1_gamma_inc_like, fmt_erfc_like, fmt1_erfc_like

   ! C: SML_FLOAT64 = DBL_EPSILON * .5
   real(rk), parameter :: sml = epsilon(1.0_rk) * 0.5_rk
   real(rk), parameter :: sqrtpie4 = &
      0.8862269254527580136490837416705725913987747280611935641069038949264_rk
   real(rk), parameter :: erfc_bound = 200.0_rk

contains

#include "cint_fmt_body.inc"

end module cint_fmt_dp

module cint_fmt_qp
   ! `rk => dp` first: LFortran 0.64 drops the plain `dp` otherwise (cint_fmt.F90).
   use cint_const,  only: rk => dp, dp
   use cint_quad, only: quad, operator(+), operator(-), operator(*), operator(/), &
                        operator(<), operator(>), operator(<=), operator(>=), &
                        operator(==), operator(**), assignment(=), quad_from, quad_to_dp, &
                        quad_from3, sqrt, abs, exp, erf, erfc
#undef RKTYPE
#undef TO_RK
#undef TO_DP
#define RKTYPE type(quad)
#define TO_RK(x) quad_from(real(x, dp))
#define TO_DP(x) quad_to_dp(x)
   use cint_fmt_tab, only: turnover_point_qp => turnover_point
   implicit none
   private
   public :: gamma_inc_like, fmt1_gamma_inc_like, fmt_erfc_like, fmt1_erfc_like

   ! C: SML_FLOAT128 = 1.0e-35.  The C also has SML_FLOAT80 = 2.0e-20 for its
   ! 80-bit ladder; that ladder becomes binary128 here, so it gets the tighter
   ! tolerance and converges further than the C does.
   real(rk), parameter :: sml = 1.0e-35_rk
   ! sqrt(pi)/2 to 67 digits -- it is the leading factor in both the Boys
   ! function and the moments, so truncating it caps the whole extended
   ! ladder.  A quad cannot be a named constant (rebuilding one is a call
   ! into the C shim), so the name is a macro and the rebuild happens at
   ! each use site: it appears once per entry call, never in an inner loop,
   ! and the module stays free of state that would need a guarded init.
   ! The three parts sum to the exact binary128 value; quad_const_check
   ! measures that, as raw bits, against gfortran's own parse of the digits.
   real(dp), parameter :: sqrtpie4_3d(3) = &
      [0.886226925452758_dp, -3.8332932499128993e-17_dp, -6.5291674539727145e-34_dp]
#define sqrtpie4 quad_from3(sqrtpie4_3d(1), sqrtpie4_3d(2), sqrtpie4_3d(3))
   real(rk), parameter :: erfc_bound = 200.0_rk

   ! The turnover table is a set of thresholds, not data feeding the result,
   ! so double precision is ample; promote on read.
   real(rk), parameter :: turnover_point(0:39) = real(turnover_point_qp, rk)

contains

#include "cint_fmt_body.inc"

end module cint_fmt_qp
