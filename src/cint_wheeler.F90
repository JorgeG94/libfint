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
   use cint_const,       only: dp, rk => dp
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
   use cint_const,          only: dp, rk => qp
   use cint_fmt_qp,         only: fmt_erfc_like
   use cint_eigh,           only: cint_diagonalize
   use cint_tab_jacobi,     only: jacobi_coef_order => JACOBI_COEF_ORDER
   use cint_tab_jacobi_ext, only: jacobi_coef     => qJACOBI_COEF, &
                                  jacobi_rn_part2 => qJACOBI_RN_PART2, &
                                  jacobi_sn       => qJACOBI_SN, &
                                  jacobi_alpha_rk => qJACOBI_ALPHA, &
                                  jacobi_beta_rk  => qJACOBI_BETA
   implicit none
   private
   public :: rys_laguerre, rys_jacobi

   integer,  parameter :: MXRYSROOTS = 32
   integer,  parameter :: flocke_extra_order = 36
   real(rk), parameter :: sqrtpie4 = &
      0.8862269254527580136490837416705725913987747280611935641069038949264_rk
   real(rk), parameter :: smallx_limit = 3.0e-7_rk

contains

#include "cint_wheeler_body.inc"

end module cint_wheeler_qp
