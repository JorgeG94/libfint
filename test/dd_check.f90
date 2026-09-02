program dd_check
   !! Does double-double actually carry ~106 bits, measured against real128?
   use iso_fortran_env, only: real64, real128
   use cint_dd, only: dd, dd_from, dd_to_dp, dd_add, dd_mul, two_sum, two_prod
   implicit none
   integer, parameter :: dp = real64, qp = real128
   integer :: i, n, seed_size
   integer, allocatable :: seed(:)
   real(dp) :: a, b, p, e, s, se
   real(qp) :: exact, got
   real(dp) :: worst_sum, worst_prod, rel
   type(dd) :: x, y, z

   call random_seed(size=seed_size)
   allocate (seed(seed_size)); seed = 20260826
   call random_seed(put=seed)

   ! 1. two_sum and two_prod must be EXACT: s + e == a + b in quad, bit for bit.
   n = 200000
   worst_sum = 0.0_dp; worst_prod = 0.0_dp
   do i = 1, n
      call random_number(a); call random_number(b)
      a = (a - 0.5_dp) * 1.0e6_dp; b = (b - 0.5_dp) * 1.0e-3_dp
      call two_sum(a, b, s, se)
      exact = real(a, qp) + real(b, qp)
      got = real(s, qp) + real(se, qp)
      if (got /= exact) worst_sum = 1.0_dp
      call two_prod(a, b, p, e)
      exact = real(a, qp) * real(b, qp)
      got = real(p, qp) + real(e, qp)
      if (got /= exact) worst_prod = 1.0_dp
   end do
   print *, "two_sum  exact on all", n, "pairs :", worst_sum == 0.0_dp
   print *, "two_prod exact on all", n, "pairs :", worst_prod == 0.0_dp
   ! Exact means exact: one inexact pair is a broken transformation, not noise.
   if (worst_sum /= 0.0_dp) error stop "two_sum is not exact"
   if (worst_prod /= 0.0_dp) error stop "two_prod is not exact"

   ! 2. Accumulated dd arithmetic against the same computation in real128.
   !    A sum of products is where a plain double loses ground fastest.
   block
      real(qp) :: acc_q
      type(dd) :: acc
      real(dp) :: acc_d, worst
      worst = 0.0_dp
      acc = dd_from(0.0_dp); acc_q = 0.0_qp; acc_d = 0.0_dp
      do i = 1, 100000
         call random_number(a); call random_number(b)
         a = a + 1.0_dp; b = b * 1.0e-8_dp
         x = dd_from(a); y = dd_from(b)
         acc = dd_add(acc, dd_mul(x, y))
         acc_q = acc_q + real(a, qp) * real(b, qp)
         acc_d = acc_d + a * b
      end do
      ! The full pair against quad -- collapsing with dd_to_dp first would only
      ! measure the final rounding to double, not what dd carried.
      rel = abs(real((real(acc%hi, qp) + real(acc%lo, qp) - acc_q) / acc_q, dp))
      print *, "dd  (hi+lo) vs quad, rel error :", rel
      print *, "dd  collapsed to double        :", &
         abs(real((real(dd_to_dp(acc), qp) - acc_q) / acc_q, dp))
      print *, "plain double vs quad, rel error:", &
         abs(real((real(acc_d, qp) - acc_q) / acc_q, dp))
      ! The point of the type. 1e-25 is far looser than the 5e-31 measured, so
      ! this fails on a broken renormalisation rather than on a rounding change.
      ! It also fails under -ffast-math, which flattens the error terms to zero
      ! and leaves double precision wearing a double-double type.
      if (rel > 1.0e-25_dp) error stop "double-double is no better than double"
   end block

   ! 3. Accuracy on a hard case: cancellation. Two nearly equal quantities
   !    subtracted is where double dies and where the Rys roots actually live.
   block
      real(qp) :: q1, q2, qref
      type(dd) :: d1, d2, ddiff
      real(dp) :: dref, ddif
      a = 1.0_dp + 2.0_dp**(-20)
      b = 1.0_dp
      d1 = dd_mul(dd_from(a), dd_from(a))
      d2 = dd_mul(dd_from(b), dd_from(b))
      ddiff = dd_add(d1, dd_mul(d2, dd_from(-1.0_dp)))
      q1 = real(a, qp)*real(a, qp); q2 = real(b, qp)*real(b, qp)
      qref = q1 - q2
      dref = a*a - b*b
      ddif = dd_to_dp(ddiff)
      print *, "a*a - b*b, a = 1+2^-20"
      print *, "   quad   :", real(qref, dp)
      print *, "   dd     : rel err", abs(real((real(ddif, qp) - qref)/qref, dp))
      print *, "   double : rel err", abs(real((real(dref, qp) - qref)/qref, dp))
   end block
end program dd_check
