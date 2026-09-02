!
! Rys roots and weights by Schmidt orthogonalisation (the RDK scheme), at both
! working precisions.
!
! As with cint_fmt, the algorithm is written once -- in cint_schmidt_body.inc
! -- and included into two modules, one bound to dp and one to qp.  Each pulls
! in the Boys function of its own precision, so the moment sequence is
! generated and orthogonalised entirely at `rk`; only the root finding drops to
! double, because _CINT_polynomial_roots has no extended variant in the C
! either.
!
! Maps onto the C as:
!
!   CINTrys_schmidt   (R_dsmit, double)       ->  cint_schmidt_dp
!   CINTlrys_schmidt  (R_lsmit, long double)  ->  cint_schmidt_qp
!   CINTqrys_schmidt  (R_qsmit, __float128)   ->  cint_schmidt_qp
!
! The C's three copies have drifted apart in two small ways, and each
! instantiation here reproduces its own counterpart rather than papering over
! the difference -- which is what makes the double path comparable to the C bit
! for bit:
!
!   * nroots == 1.  The double copy returns f1/(f0-f1) directly; the extended
!     copies form f1/f0 and let the shared tail apply root/(1-root).  Equal on
!     paper, but they round differently, at the 1e-15 level.
!   * r_smit failure.  The double copy normalises the code to 1; the extended
!     copies return the order at which orthogonality was lost.  No caller looks
!     at anything but zero/non-zero.
!
module cint_schmidt_dp
   ! `rk => dp` first: LFortran 0.64 drops the plain `dp` otherwise (cint_fmt.F90).
   use cint_const,      only: rk => dp, dp
#undef RKTYPE
#undef TO_RK
#undef TO_DP
#define RKTYPE real(rk)
#define TO_RK(x) real(x, rk)
#define TO_DP(x) real(x, dp)
   use cint_fmt_dp,     only: gamma_inc_like, fmt_erfc_like
   use cint_find_roots, only: cint_polynomial_roots, poly_value1
   implicit none
   private
   public :: rys_schmidt, r_smit

   integer, parameter :: MXRYSROOTS = 32
   logical, parameter :: one_root_early_return = .true.
   logical, parameter :: smit_err_as_one       = .true.

contains

#include "cint_schmidt_body.inc"

end module cint_schmidt_dp

module cint_schmidt_qp
   ! `rk => dp` first: LFortran 0.64 drops the plain `dp` otherwise (cint_fmt.F90).
   use cint_const,      only: rk => dp, dp
   use cint_quad, only: quad, operator(+), operator(-), operator(*), operator(/), &
                        operator(<), operator(>), operator(<=), operator(>=), &
                        operator(==), operator(**), assignment(=), quad_from, quad_to_dp, &
                        sqrt, abs, exp, erf, erfc
#undef RKTYPE
#undef TO_RK
#undef TO_DP
#define RKTYPE type(quad)
#define TO_RK(x) quad_from(real(x, dp))
#define TO_DP(x) quad_to_dp(x)
   use cint_fmt_qp,     only: gamma_inc_like, fmt_erfc_like
   use cint_find_roots, only: cint_polynomial_roots, poly_value1
   implicit none
   private
   public :: rys_schmidt, r_smit

   integer, parameter :: MXRYSROOTS = 32
   logical, parameter :: one_root_early_return = .false.
   logical, parameter :: smit_err_as_one       = .false.

contains

#include "cint_schmidt_body.inc"

end module cint_schmidt_qp
