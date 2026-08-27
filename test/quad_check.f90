program q_check
   !! The shim must reproduce real128 EXACTLY -- gfortran's real128 is
   !! libquadmath, and so is this. Anything short of bit-equality means the
   !! plumbing is lossy, because the arithmetic underneath is the same code.
   !!
   !! Compared as raw bits, all 128 of them: converting to double first would
   !! forgive 60 low bits of drift, which is precisely the kind of loss this
   !! test exists to catch.
   use iso_fortran_env, only: real64, real128, int64
   use cint_quad, only: quad, quad_from, operator(+), operator(-), &
                        operator(*), operator(/), sqrt, exp, erf, erfc
   implicit none
   integer, parameter :: dp = real64, qp = real128
   integer :: i, nbad
   real(dp) :: x
   real(qp) :: ref
   type(quad) :: z

   nbad = 0
   do i = 1, 2000
      x = real(i, dp) * 0.0031_dp
      z = sqrt(quad_from(x));  ref = sqrt(real(x, qp));   if (differs(z, ref)) nbad = nbad + 1
      z = exp(quad_from(x));   ref = exp(real(x, qp));    if (differs(z, ref)) nbad = nbad + 1
      z = erf(quad_from(x));   ref = erf(real(x, qp));    if (differs(z, ref)) nbad = nbad + 1
      z = erfc(quad_from(x));  ref = erfc(real(x, qp));   if (differs(z, ref)) nbad = nbad + 1
      z = quad_from(x)*quad_from(x) - quad_from(1.0_dp)
      ref = real(x, qp)*real(x, qp) - 1.0_qp
      if (differs(z, ref)) nbad = nbad + 1
   end do
   print '(a,i0,a)', "bitwise mismatches vs real128 over 10000 comparisons: ", nbad, ""
   if (nbad /= 0) error stop "C shim does not reproduce real128"
   print *, "OK -- bit-identical to real128"

contains

   logical function differs(a, b) result(v)
      type(quad), intent(in) :: a
      real(qp), intent(in) :: b
      v = any(transfer(a, [0_int64, 0_int64]) /= transfer(b, [0_int64, 0_int64]))
   end function differs

end program q_check
