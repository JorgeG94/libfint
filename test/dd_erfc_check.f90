program dd_erfc_t
   use iso_fortran_env, only: real64, real128
   use cint_dd, only: dd, dd_from, erfc, operator(-)
   implicit none
   integer, parameter :: dp = real64, qp = real128
   integer :: i
   real(dp) :: x, rel, worst(5)
   real(qp) :: qref
   type(dd) :: z
   character(len=18) :: names(5)
   names = [character(len=18) :: "A  x<=0.46875", "B  0.46875..2", "C  2..9", &
                                 "D  9..26.6", "negative x"]
   worst = 0.0_dp
   do i = 1, 5600
      x = real(i, dp) * 0.005_dp
      z = erfc(dd_from(x)); qref = erfc(real(x, qp))
      ! **Double-double degrades to double near the underflow floor, and it
      ! is the LOW word that runs out first.** A result of 1e-300 needs `lo`
      ! at ~1e-316, which is denormal and carries almost no bits. Full 106-bit
      ! precision therefore stops at hi > tiny*2**53 ~ 2e-292, not at tiny.
      ! Below that the answer is still right to double precision; below tiny
      ! it is zero. Both are checked, with the budget each can actually meet.
      !
      ! This band needs denormals to be live in the process, which is not the
      ! default everywhere -- see the flag on this target in test/CMakeLists.txt.
      if (qref <= real(tiny(1.0_dp), qp)) cycle
      if (qref < real(tiny(1.0_dp), qp) * 2.0_qp**53) then
         rel = abs(real((real(z%hi, qp) + real(z%lo, qp) - qref)/qref, dp))
         if (rel > 1.0e-15_dp) error stop "dd_erfc worse than double near underflow"
         cycle
      end if
      rel = abs(real((real(z%hi, qp) + real(z%lo, qp) - qref)/qref, dp))
      if (x <= 0.46875_dp) then;     worst(1) = max(worst(1), rel)
      else if (x <= 2.0_dp) then;    worst(2) = max(worst(2), rel)
      else if (x <= 9.0_dp) then;    worst(3) = max(worst(3), rel)
      else;                          worst(4) = max(worst(4), rel)
      end if
   end do
   do i = 1, 400
      x = -real(i, dp) * 0.01_dp
      z = erfc(dd_from(x)); qref = erfc(real(x, qp))
      rel = abs(real((real(z%hi, qp) + real(z%lo, qp) - qref)/qref, dp))
      worst(5) = max(worst(5), rel)
   end do
   do i = 1, 5
      print '(a,a,es12.4)', "  ", names(i), worst(i)
   end do
   if (maxval(worst) > 1.0e-28_dp) error stop "dd_erfc is not accurate enough"
   print *, "OK"
end program dd_erfc_t
