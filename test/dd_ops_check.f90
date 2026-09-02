program dd_ops
   !! div and sqrt against real128, plus the operator-overload path the
   !! shared bodies will actually take.
   use iso_fortran_env, only: real64, real128
   use cint_dd, only: dd, dd_from, dd_to_dp, operator(+), operator(-), &
                      operator(*), operator(/), operator(<), sqrt, abs
   implicit none
   integer, parameter :: dp = real64, qp = real128
   integer :: i
   real(dp) :: a, b, worst_div, worst_sqrt, rel
   real(qp) :: qa, qb, qref
   type(dd) :: x, y, z

   worst_div = 0.0_dp; worst_sqrt = 0.0_dp
   do i = 1, 50000
      call random_number(a); call random_number(b)
      a = a * 1.0e3_dp + 1.0e-3_dp
      b = b * 1.0e3_dp + 1.0e-3_dp
      x = dd_from(a); y = dd_from(b)

      z = x / y
      qa = real(a, qp); qb = real(b, qp); qref = qa / qb
      rel = abs(real((real(z%hi, qp) + real(z%lo, qp) - qref) / qref, dp))
      worst_div = max(worst_div, rel)

      z = sqrt(x)
      qref = sqrt(qa)
      rel = abs(real((real(z%hi, qp) + real(z%lo, qp) - qref) / qref, dp))
      worst_sqrt = max(worst_sqrt, rel)
   end do
   print *, "dd_div  worst relative error :", worst_div
   print *, "dd_sqrt worst relative error :", worst_sqrt

   ! The expression form the .inc bodies use, mixing dd and literals.
   block
      type(dd) :: tt, f0
      real(qp) :: qtt, qf0
      tt = sqrt(dd_from(2.0_dp))
      f0 = 0.5_dp * tt / (tt + 1.0_dp) - tt * 0.25_dp
      qtt = sqrt(2.0_qp)
      qf0 = 0.5_qp * qtt / (qtt + 1.0_qp) - qtt * 0.25_qp
      rel = abs(real((real(f0%hi, qp) + real(f0%lo, qp) - qf0) / qf0, dp))
      print *, "mixed dd/literal expression  :", rel
      if (rel > 1.0e-28_dp) error stop "operator overloads lose precision"
   end block

   if (worst_div > 1.0e-28_dp) error stop "dd_div is not accurate enough"
   if (worst_sqrt > 1.0e-28_dp) error stop "dd_sqrt is not accurate enough"
   print *, "OK"
end program dd_ops
