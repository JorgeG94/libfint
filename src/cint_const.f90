!
! Kind parameters for the Fortran port.
!
module cint_const
   use iso_fortran_env, only: real64, real128
   implicit none
   public
   integer, parameter :: dp = real64
   ! The extended ladder.  PORT_TO_FORTRAN.md 3.1: the C runs two ladders,
   ! 80-bit long double and __float128; both collapse onto binary128 here,
   ! which is simpler and strictly more accurate.  gfortran and ifx have it;
   ! nvfortran does not, which is the constraint to remember.
   integer, parameter :: qp = real128
end module cint_const
