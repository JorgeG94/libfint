program dd_fn
   !! exp and erf in double-double, against real128 in the same binary.
   use iso_fortran_env, only: real64, real128
   use cint_dd, only: dd, dd_from, operator(-), exp, erf
   implicit none
   integer, parameter :: dp = real64, qp = real128
   integer :: i
   real(dp) :: x, worst_exp, worst_erf, rel
   real(qp) :: qref
   type(dd) :: z

   worst_exp = 0.0_dp
   do i = 0, 1400
      x = -35.0_dp + real(i, dp) * 0.05_dp
      z = exp(dd_from(x))
      qref = exp(real(x, qp))
      if (qref /= 0.0_qp) then
         rel = abs(real((real(z%hi, qp) + real(z%lo, qp) - qref) / qref, dp))
         worst_exp = max(worst_exp, rel)
      end if
   end do
   print *, "dd_exp worst rel error, x in [-35,35] :", worst_exp

   ! erf over the range the Boys function actually asks for.
   worst_erf = 0.0_dp
   do i = 0, 1300
      x = real(i, dp) * 0.005_dp
      z = erf(dd_from(x))
      qref = erf(real(x, qp))
      if (qref /= 0.0_qp) then
         rel = abs(real((real(z%hi, qp) + real(z%lo, qp) - qref) / qref, dp))
         worst_erf = max(worst_erf, rel)
      end if
   end do
   print *, "dd_erf worst rel error, x in [0,6.5]  :", worst_erf

   ! 1e-28 against 3.3e-30 measured: fails on a broken reduction or a bad
   ! constant, not on a rounding change. A wrong DD_LN2 lands at 1.4e-17.
   if (worst_exp > 1.0e-28_dp) error stop "dd_exp is not accurate enough"
   if (worst_erf > 1.0e-28_dp) error stop "dd_erf is not accurate enough"
   print *, "OK"
end program dd_fn
