program t3
   use iso_fortran_env, only: real64, real128
   use cint_quad, only: quad, quad_from3, quad_to_dp, operator(-), operator(/), abs
   implicit none
   integer, parameter :: dp = real64, qp = real128
   type(quad) :: q
   real(qp) :: ref
   ! sqrt(pi)/2 to 67 digits, split three ways
   q = quad_from3(0.886226925452758_dp, -3.8332932499128993e-17_dp, -6.5291674539727145e-34_dp)
   ref = 0.8862269254527580136490837416705725913987747280611935641069038949264_qp
   print '(a,es24.16)', "reconstructed - real128 (as double) = ", quad_to_dp(q) - real(ref, dp)
   if (quad_to_dp(q) /= real(ref, dp)) error stop "three-double split is not exact"
   print *, "OK -- 67-digit constant reaches quad exactly"
end program
