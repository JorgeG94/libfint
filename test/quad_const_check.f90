program t3
   !! quad_from3 must land a 67-digit constant on the SAME binary128 that
   !! gfortran's own literal parsing produces -- checked as raw bits.  Two
   !! doubles are seven bits short of the 113-bit significand; three are
   !! exact, and this is where that claim is measured rather than asserted.
   use iso_fortran_env, only: real64, real128, int64
   use cint_quad, only: quad, quad_from3, quad_to_dp
   implicit none
   integer, parameter :: dp = real64, qp = real128
   type(quad) :: q
   real(qp) :: ref
   ! sqrt(pi)/2 to 67 digits, split three ways against exact binary remainders
   q = quad_from3(0.886226925452758_dp, -3.8332932499128993e-17_dp, -6.5291674539727145e-34_dp)
   ref = 0.8862269254527580136490837416705725913987747280611935641069038949264_qp
   print '(a,es24.16)', "reconstructed - real128 (as double) = ", quad_to_dp(q) - real(ref, dp)
   if (any(transfer(q, [0_int64, 0_int64]) /= transfer(ref, [0_int64, 0_int64]))) then
      error stop "three-double split does not hit the real128 literal"
   end if
   print *, "OK -- 67-digit constant reaches quad bit-exactly"
end program
