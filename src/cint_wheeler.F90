!
! Rys roots and weights by the Wheeler recursion, at both working precisions.
!
! Same include-once shape as cint_fmt and cint_schmidt.  Each instantiation
! binds the moment tables of its own precision:
!
!   cint_wheeler_dp  ->  JACOBI_*    (double)
!   cint_wheeler_qp  ->  qJACOBI_*   (binary128)
!
! The C has a third set, lJACOBI_* at 80-bit, which collapses into the
! quadruple instantiation here along with everything else on the extended
! ladder (PORT_TO_FORTRAN.md 3.1).  The Flocke seed order follows suit: the C
! uses 20 extra orders for double, 24 for long double and 36 for __float128,
! and the quadruple instantiation takes 36.
!
module cint_wheeler_dp
   ! `rk => dp` first: LFortran 0.64 drops the plain `dp` otherwise (cint_fmt.F90).
   use cint_const,       only: rk => dp, dp
#undef RKTYPE
#undef TO_RK
#undef TO_DP
#undef ENSURE_EXT_TABLES
#define ENSURE_EXT_TABLES
#define RKTYPE real(rk)
#define TO_RK(x) real(x, rk)
#define TO_DP(x) real(x, dp)
   use cint_fmt_dp,      only: fmt_erfc_like
   use cint_eigh,        only: cint_diagonalize
   use cint_tab_jacobi,  only: jacobi_coef       => JACOBI_COEF, &
                               jacobi_coef_order => JACOBI_COEF_ORDER, &
                               jacobi_rn_part2   => JACOBI_RN_PART2, &
                               jacobi_sn         => JACOBI_SN, &
                               jacobi_alpha_rk   => JACOBI_ALPHA, &
                               jacobi_beta_rk    => JACOBI_BETA
   implicit none
   private
   public :: rys_laguerre, rys_jacobi

   integer,  parameter :: MXRYSROOTS = 32
   integer,  parameter :: flocke_extra_order = 20
   real(rk), parameter :: sqrtpie4 = &
      0.8862269254527580136490837416705725913987747280611935641069038949264_rk
   real(rk), parameter :: smallx_limit = 3.0e-7_rk

contains

#include "cint_wheeler_body.inc"

end module cint_wheeler_dp

module cint_wheeler_qp
   ! `rk => dp` first: LFortran 0.64 drops the plain `dp` otherwise (cint_fmt.F90).
   use cint_const,          only: rk => dp, dp
   use cint_quad, only: quad, operator(+), operator(-), operator(*), operator(/), &
                        operator(<), operator(>), operator(<=), operator(>=), &
                        operator(==), operator(**), assignment(=), quad_from, quad_to_dp, &
                        quad_from3, sqrt, abs, exp, erf, erfc
#undef RKTYPE
#undef TO_RK
#undef TO_DP
#undef ENSURE_EXT_TABLES
#define RKTYPE type(quad)
#define TO_RK(x) quad_from(real(x, dp))
#define TO_DP(x) quad_to_dp(x)
! The quad tables cannot be parameters (see cint_tab_jacobi_ext), so the one
! routine that reads them has to fill them first; the double instantiation
! defines this hook empty because its tables are compile-time constants.
#define ENSURE_EXT_TABLES call cint_tab_jacobi_ext_init()
   use cint_fmt_qp,         only: fmt_erfc_like
   use cint_eigh,           only: cint_diagonalize
   use cint_tab_jacobi,     only: jacobi_coef_order => JACOBI_COEF_ORDER
   use cint_tab_jacobi_ext, only: cint_tab_jacobi_ext_init, &
                                  jacobi_coef     => qJACOBI_COEF, &
                                  jacobi_rn_part2 => qJACOBI_RN_PART2, &
                                  jacobi_sn       => qJACOBI_SN, &
                                  jacobi_alpha_rk => qJACOBI_ALPHA, &
                                  jacobi_beta_rk  => qJACOBI_BETA
   implicit none
   private
   public :: rys_laguerre, rys_jacobi

   integer,  parameter :: MXRYSROOTS = 32
   integer,  parameter :: flocke_extra_order = 36
   ! sqrt(pi)/2 to 67 digits -- the leading factor in the Boys function and
   ! the moments, so truncating it caps the whole extended ladder.  A quad
   ! cannot be a named constant, so the name is a macro rebuilding it at
   ! each use site through quad_from3; once per entry call, never in a loop.
   real(dp), parameter :: sqrtpie4_3d(3) = &
      [0.886226925452758_dp, -3.8332932499128993e-17_dp, -6.5291674539727145e-34_dp]
#define sqrtpie4 quad_from3(sqrtpie4_3d(1), sqrtpie4_3d(2), sqrtpie4_3d(3))
   real(rk), parameter :: smallx_limit = 3.0e-7_rk

contains

#include "cint_wheeler_body.inc"

end module cint_wheeler_qp
