!
! binary128 through C, for compilers with no real128.
!
! **This mirrors cint_dd's interface exactly** -- same operators, same generic
! names, same conversion helpers -- so the extended instantiation selects
! between them by changing which module it uses, and nothing else. That is the
! seam the RKTYPE macro created.
!
! Where cint_dd trades 7 bits of significand for being all-Fortran and fast,
! this keeps all 113 and reproduces what gfortran's real128 produces today,
! because gfortran's real128 *is* libquadmath. It costs a call per operation.
!
module cint_quad
   use iso_c_binding, only: c_double, c_int
   use cint_const, only: dp
   implicit none
   private

   type, bind(c), public :: quad
      !! 16 bytes, aligned as a double pair. Never inspected on the Fortran
      !! side -- only the C shim knows what the bits mean.
      !!
      !! `bind(c)` is required for the interfaces below to be interoperable,
      !! and it forbids default initialisation -- so unlike cint_dd's type
      !! this one starts undefined and every use must assign first.
      real(c_double) :: w(2)
   end type quad

   public :: quad_from, quad_to_dp
   public :: operator(+), operator(-), operator(*), operator(/)
   public :: operator(<), operator(>), operator(<=), operator(>=), operator(==)
   public :: sqrt, abs, exp, erf, erfc, max, min, assignment(=), operator(**)

   interface quad_from
      module procedure q_from_r, q_from_i
   end interface
   interface operator(+)
      module procedure q_add, q_add_r, q_r_add
   end interface
   interface operator(-)
      module procedure q_sub, q_sub_r, q_r_sub, q_neg
   end interface
   interface operator(*)
      module procedure q_mul, q_mul_r, q_r_mul
   end interface
   interface operator(/)
      module procedure q_div, q_div_r, q_r_div
   end interface
   interface operator(<)
      module procedure q_lt, q_lt_r
   end interface
   interface operator(>)
      module procedure q_gt, q_gt_r
   end interface
   interface operator(<=)
      module procedure q_le, q_le_r
   end interface
   interface operator(>=)
      module procedure q_ge, q_ge_r
   end interface
   interface operator(==)
      module procedure q_eq, q_eq_r
   end interface
   interface operator(**)
      module procedure q_pow_i
   end interface
   interface assignment(=)
      module procedure q_assign_r
   end interface
   interface sqrt
      module procedure q_sqrt
   end interface
   interface abs
      module procedure q_abs
   end interface
   interface exp
      module procedure q_exp
   end interface
   interface erf
      module procedure q_erf
   end interface
   interface erfc
      module procedure q_erfc
   end interface
   interface max
      module procedure q_max, q_max_r
   end interface
   interface min
      module procedure q_min, q_min_r
   end interface

   interface
      subroutine mqq_add(a, b, r) bind(c, name="mqq_add")
         import quad
         type(quad), intent(in) :: a, b
         type(quad), intent(out) :: r
      end subroutine mqq_add
      subroutine mqq_sub(a, b, r) bind(c, name="mqq_sub")
         import quad
         type(quad), intent(in) :: a, b
         type(quad), intent(out) :: r
      end subroutine mqq_sub
      subroutine mqq_mul(a, b, r) bind(c, name="mqq_mul")
         import quad
         type(quad), intent(in) :: a, b
         type(quad), intent(out) :: r
      end subroutine mqq_mul
      subroutine mqq_div(a, b, r) bind(c, name="mqq_div")
         import quad
         type(quad), intent(in) :: a, b
         type(quad), intent(out) :: r
      end subroutine mqq_div
      subroutine mqq_neg(a, r) bind(c, name="mqq_neg")
         import quad
         type(quad), intent(in) :: a
         type(quad), intent(out) :: r
      end subroutine mqq_neg
      subroutine mqq_abs(a, r) bind(c, name="mqq_abs")
         import quad
         type(quad), intent(in) :: a
         type(quad), intent(out) :: r
      end subroutine mqq_abs
      subroutine mqq_sqrt(a, r) bind(c, name="mqq_sqrt")
         import quad
         type(quad), intent(in) :: a
         type(quad), intent(out) :: r
      end subroutine mqq_sqrt
      subroutine mqq_exp(a, r) bind(c, name="mqq_exp")
         import quad
         type(quad), intent(in) :: a
         type(quad), intent(out) :: r
      end subroutine mqq_exp
      subroutine mqq_erf(a, r) bind(c, name="mqq_erf")
         import quad
         type(quad), intent(in) :: a
         type(quad), intent(out) :: r
      end subroutine mqq_erf
      subroutine mqq_erfc(a, r) bind(c, name="mqq_erfc")
         import quad
         type(quad), intent(in) :: a
         type(quad), intent(out) :: r
      end subroutine mqq_erfc
      subroutine mqq_from_d(d, r) bind(c, name="mqq_from_d")
         import quad, c_double
         real(c_double), value :: d
         type(quad), intent(out) :: r
      end subroutine mqq_from_d
      subroutine mqq_from_i(n, r) bind(c, name="mqq_from_i")
         import quad, c_int
         integer(c_int), value :: n
         type(quad), intent(out) :: r
      end subroutine mqq_from_i
      function mqq_to_d(a) bind(c, name="mqq_to_d") result(v)
         import quad, c_double
         type(quad), intent(in) :: a
         real(c_double) :: v
      end function mqq_to_d
      function mqq_cmp(a, b) bind(c, name="mqq_cmp") result(v)
         import quad, c_int
         type(quad), intent(in) :: a, b
         integer(c_int) :: v
      end function mqq_cmp
   end interface

contains

   type(quad) function q_add(a, b) result(r)
      type(quad), intent(in) :: a, b
      call mqq_add(a, b, r)
   end function q_add

   type(quad) function q_sub(a, b) result(r)
      type(quad), intent(in) :: a, b
      call mqq_sub(a, b, r)
   end function q_sub

   type(quad) function q_mul(a, b) result(r)
      type(quad), intent(in) :: a, b
      call mqq_mul(a, b, r)
   end function q_mul

   type(quad) function q_div(a, b) result(r)
      type(quad), intent(in) :: a, b
      call mqq_div(a, b, r)
   end function q_div

   type(quad) function q_neg(a) result(r)
      type(quad), intent(in) :: a
      call mqq_neg(a, r)
   end function q_neg

   type(quad) function q_abs(a) result(r)
      type(quad), intent(in) :: a
      call mqq_abs(a, r)
   end function q_abs

   type(quad) function q_sqrt(a) result(r)
      type(quad), intent(in) :: a
      call mqq_sqrt(a, r)
   end function q_sqrt

   type(quad) function q_exp(a) result(r)
      type(quad), intent(in) :: a
      call mqq_exp(a, r)
   end function q_exp

   type(quad) function q_erf(a) result(r)
      type(quad), intent(in) :: a
      call mqq_erf(a, r)
   end function q_erf

   type(quad) function q_erfc(a) result(r)
      type(quad), intent(in) :: a
      call mqq_erfc(a, r)
   end function q_erfc

   type(quad) function q_from_r(x) result(r)
      real(dp), intent(in) :: x
      call mqq_from_d(real(x, c_double), r)
   end function q_from_r

   type(quad) function q_from_i(n) result(r)
      integer, intent(in) :: n
      call mqq_from_i(int(n, c_int), r)
   end function q_from_i

   real(dp) function quad_to_dp(a) result(v)
      type(quad), intent(in) :: a
      v = real(mqq_to_d(a), dp)
   end function quad_to_dp

   subroutine q_assign_r(out, in)
      type(quad), intent(out) :: out
      real(dp), intent(in) :: in
      out = q_from_r(in)
   end subroutine q_assign_r

   type(quad) function q_pow_i(a, n) result(r)
      type(quad), intent(in) :: a
      integer, intent(in) :: n
      type(quad) :: base
      integer :: k
      r = q_from_r(1.0_dp)
      if (n == 0) return
      base = a; k = abs(n)
      do while (k > 0)
         if (mod(k, 2) == 1) r = q_mul(r, base)
         base = q_mul(base, base)
         k = k/2
      end do
      if (n < 0) r = q_div(q_from_r(1.0_dp), r)
   end function q_pow_i

   type(quad) function q_max(a, b) result(r)
      type(quad), intent(in) :: a, b
      r = a
      if (mqq_cmp(b, a) > 0) r = b
   end function q_max

   type(quad) function q_min(a, b) result(r)
      type(quad), intent(in) :: a, b
      r = a
      if (mqq_cmp(b, a) < 0) r = b
   end function q_min

   type(quad) function q_add_r(a, b) result(r)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      r = q_add(a, q_from_r(b))
   end function q_add_r

   type(quad) function q_r_add(a, b) result(r)
      real(dp), intent(in) :: a
      type(quad), intent(in) :: b
      r = q_add(q_from_r(a), b)
   end function q_r_add

   type(quad) function q_sub_r(a, b) result(r)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      r = q_sub(a, q_from_r(b))
   end function q_sub_r

   type(quad) function q_r_sub(a, b) result(r)
      real(dp), intent(in) :: a
      type(quad), intent(in) :: b
      r = q_sub(q_from_r(a), b)
   end function q_r_sub

   type(quad) function q_mul_r(a, b) result(r)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      r = q_mul(a, q_from_r(b))
   end function q_mul_r

   type(quad) function q_r_mul(a, b) result(r)
      real(dp), intent(in) :: a
      type(quad), intent(in) :: b
      r = q_mul(q_from_r(a), b)
   end function q_r_mul

   type(quad) function q_div_r(a, b) result(r)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      r = q_div(a, q_from_r(b))
   end function q_div_r

   type(quad) function q_r_div(a, b) result(r)
      real(dp), intent(in) :: a
      type(quad), intent(in) :: b
      r = q_div(q_from_r(a), b)
   end function q_r_div

   logical function q_lt(a, b) result(v)
      type(quad), intent(in) :: a, b
      v = mqq_cmp(a, b) < 0
   end function q_lt

   logical function q_lt_r(a, b) result(v)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      v = q_lt(a, q_from_r(b))
   end function q_lt_r

   logical function q_gt(a, b) result(v)
      type(quad), intent(in) :: a, b
      v = mqq_cmp(a, b) > 0
   end function q_gt

   logical function q_gt_r(a, b) result(v)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      v = q_gt(a, q_from_r(b))
   end function q_gt_r

   logical function q_le(a, b) result(v)
      type(quad), intent(in) :: a, b
      v = mqq_cmp(a, b) <= 0
   end function q_le

   logical function q_le_r(a, b) result(v)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      v = q_le(a, q_from_r(b))
   end function q_le_r

   logical function q_ge(a, b) result(v)
      type(quad), intent(in) :: a, b
      v = mqq_cmp(a, b) >= 0
   end function q_ge

   logical function q_ge_r(a, b) result(v)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      v = q_ge(a, q_from_r(b))
   end function q_ge_r

   logical function q_eq(a, b) result(v)
      type(quad), intent(in) :: a, b
      v = mqq_cmp(a, b) == 0
   end function q_eq

   logical function q_eq_r(a, b) result(v)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      v = q_eq(a, q_from_r(b))
   end function q_eq_r

   type(quad) function q_max_r(a, b) result(r)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      r = q_max(a, q_from_r(b))
   end function q_max_r

   type(quad) function q_min_r(a, b) result(r)
      type(quad), intent(in) :: a
      real(dp), intent(in) :: b
      r = q_min(a, q_from_r(b))
   end function q_min_r

end module cint_quad
