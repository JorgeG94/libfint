program crossbit_dd
   !! Raw-bit dump of dd results, for cross-compiler comparison.
   !!
   !! Bits, not decimals: printed values hide denormals on nvfortran and round
   !! away the `lo` term everywhere. (LANDMINES.md #12.)
   use iso_fortran_env, only: real64, int64
   use cint_dd, only: dd, dd_from, operator(+), operator(-), operator(*), &
                      operator(/), sqrt, exp, erf
   implicit none
   integer, parameter :: dp = real64
   integer :: i
   real(dp) :: x
   type(dd) :: z
   integer(int64) :: fold_hi, fold_lo

   fold_hi = 0_int64; fold_lo = 0_int64
   do i = 1, 400
      x = real(i, dp) * 0.017_dp
      z = sqrt(dd_from(x))
      call fold(z)
      z = exp(dd_from(x - 3.0_dp))
      call fold(z)
      z = erf(dd_from(x))
      call fold(z)
      z = dd_from(1.0_dp) / dd_from(x)
      call fold(z)
      z = dd_from(x) * dd_from(x) - dd_from(1.0_dp)
      call fold(z)
   end do
   print '(a,z16.16)', "fold(hi) = 0x", fold_hi
   print '(a,z16.16)', "fold(lo) = 0x", fold_lo

contains
   subroutine fold(v)
      type(dd), intent(in) :: v
      ! Salted rotate: an unsalted xor-fold is blind to repeated values.
      fold_hi = ieor(ishftc(fold_hi, 7), transfer(v%hi, 1_int64))
      fold_lo = ieor(ishftc(fold_lo, 11), transfer(v%lo, 1_int64))
   end subroutine fold
end program crossbit_dd
