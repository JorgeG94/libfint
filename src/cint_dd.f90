!
! Double-double arithmetic: ~106 bits of significand out of two real64s.
!
! Why this exists.  The extended ladder in cint_const uses real128, and
! nvfortran and LLVM Flang do not have it -- `iso_fortran_env::real128` is -1
! there and `real(real128)` is a compile error, so libfint cannot be built with
! either.  The quad precision itself is not the problem: nvc's __float128 is a
! true binary128 with a 113-bit significand, measured, identical to gcc's.  The
! gap is the Fortran front end alone.
!
! A C shim cannot fix that, because what is missing is a *kind parameter*: every
! declaration, literal and intrinsic would have to become a call.  Two real64s
! carrying a value and its rounding error do fix it, in plain Fortran that every
! compiler has.
!
! **The algorithms are exact, and that is the point.**  `two_sum` and
! `two_prod` return the rounded result *and* the error it discarded, with no
! approximation: the error of a floating-point sum is itself representable.
! Dekker splitting is used for the product rather than an FMA, because Fortran
! has no portable FMA before F2018's IEEE_FMA and nvfortran's support for it is
! not something to rely on here.
!
! **The flag contract, measured rather than assumed.**  Folding the raw bits of
! 2000 dd results over sqrt/exp/erf/div/mul:
!
!   gfortran -O2                     hi 0854955B40D86FC2   lo 85F4EDB9AB7A5253
!   nvfortran -O2 -Mnofma            hi 0854955B40D86FC2   lo 85F4EDB9AB7A5253
!   nvfortran -O2   (FMA is default) hi 0854955B40D86FC2   lo B7E34D1CF97CFA3B
!   nvfortran -fast                  hi 0854955B40D86FC2   lo B7E34D1CF97CFA3B
!   gfortran -O2 -ffast-math         hi 1B7A0DE5351EEF18   lo 0000000000000000
!
! Two different failures, and they need different responses.
!
! **FMA contraction costs reproducibility, not accuracy.**  `hi` is identical
! everywhere; only the error term changes, and it stays a valid error term --
! accuracy over [-35,35] is 3.2e-30 contracted against 3.3e-30 not.  So an FMA
! build is *correct* and merely disagrees bit-for-bit with every other build.
! nvfortran contracts by default at -O2, so -Mnofma is required to match.  On
! ifx the spelling is -no-fma; -fno-fma is accepted and silently ignored.
!
! **-ffast-math destroys it outright.**  Every `lo` folds to exactly zero: the
! reassociation of `(a + b) - a` that the error terms are built from is gone,
! and what is left is double precision wearing a double-double type.  That is
! the one that must never be allowed, and dd_check catches it.
!
! Both findings are from this repository's own compilers, not from lore; the
! landmine catalogue in bitwise_adventures documents the same FMA behaviour
! from an independent direction.
!
module cint_dd
   use iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter :: dp = real64

   type, public :: dd
      real(dp) :: hi = 0.0_dp
      real(dp) :: lo = 0.0_dp
   end type dd

   public :: dd_from, dd_to_dp, two_sum, two_prod, dd_add, dd_mul
   public :: dd_sub, dd_div, dd_sqrt
   public :: operator(+), operator(-), operator(*), operator(/)
   public :: operator(<), operator(>), operator(<=), operator(>=), operator(==)
   public :: sqrt, abs
   public :: dd_exp, dd_erf
   public :: exp, erf

   interface exp
      module procedure dd_exp
   end interface
   interface erf
      module procedure dd_erf
   end interface

   ! Split so `hi + lo` reproduces the constant to ~1e-33. Hard-coded rather
   ! than derived from real128 at startup, because real128 is precisely what
   ! this module exists to avoid needing.
   !
   ! **`lo` is the remainder against the exact binary value of `hi`, not against
   ! its printed form.** Generating these from `repr(hi)` -- the shortest decimal
   ! that round-trips -- leaves `lo` wrong by the gap between the two, about
   ! 2e-17. The first version did exactly that, and the symptom was oddly
   ! specific: exp was accurate to 1e-29 while the range reduction chose k = 0
   ! and fell to 1.4e-17 the moment it chose anything else, because that is when
   ! LN2 first enters the arithmetic.
   type(dd), parameter :: DD_LN2 = dd(0.6931471805599453_dp, 2.3190468138462996e-17_dp)
   type(dd), parameter :: DD_INV_LN2 = dd(1.4426950408889634_dp, 2.0355273740931033e-17_dp)
   type(dd), parameter :: DD_2_SQRTPI = dd(1.1283791670955126_dp, 1.533545961316588e-17_dp)

   ! **Overloaded so the shared bodies do not change.** cint_fmt_body.inc and
   ! cint_wheeler_body.inc are included once per precision with `rk` bound to a
   ! kind; binding a derived type instead only works if `a*b - c/d` still reads
   ! as arithmetic. Everything below exists so that the algorithm text stays
   ! one copy rather than becoming two that must be kept in step.
   interface operator(+)
      module procedure dd_add, dd_add_r, dd_r_add
   end interface
   interface operator(-)
      module procedure dd_sub, dd_sub_r, dd_r_sub, dd_neg
   end interface
   interface operator(*)
      module procedure dd_mul, dd_mul_r, dd_r_mul
   end interface
   interface operator(/)
      module procedure dd_div, dd_div_r, dd_r_div
   end interface
   interface operator(<)
      module procedure dd_lt, dd_lt_r
   end interface
   interface operator(>)
      module procedure dd_gt, dd_gt_r
   end interface
   interface operator(<=)
      module procedure dd_le, dd_le_r
   end interface
   interface operator(>=)
      module procedure dd_ge, dd_ge_r
   end interface
   interface operator(==)
      module procedure dd_eq, dd_eq_r
   end interface
   interface sqrt
      module procedure dd_sqrt
   end interface
   interface abs
      module procedure dd_abs
   end interface

   ! Dekker's splitting constant, 2**27 + 1: splits a 53-bit significand into
   ! two 26-bit halves whose product is exact.
   real(dp), parameter :: SPLITTER = 134217729.0_dp

contains

   pure type(dd) function dd_from(x) result(r)
      !! A double promoted exactly: the error term is zero by construction.
      real(dp), intent(in) :: x
      r%hi = x
      r%lo = 0.0_dp
   end function dd_from

   pure real(dp) function dd_to_dp(a) result(r)
      !! Back to double, correctly rounded by the representation itself.
      type(dd), intent(in) :: a
      r = a%hi + a%lo
   end function dd_to_dp

   pure subroutine two_sum(a, b, s, e)
      !! s = fl(a+b) and e = the exact discarded error, so s + e == a + b.
      !!
      !! Knuth's version: no assumption about which operand is larger, at the
      !! cost of three more flops than Dekker's fast_two_sum.
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: s, e
      real(dp) :: bb

      s = a + b
      bb = s - a
      e = (a - (s - bb)) + (b - bb)
   end subroutine two_sum

   pure subroutine split(a, hi, lo)
      !! a = hi + lo exactly, each with a 26-bit significand.
      real(dp), intent(in) :: a
      real(dp), intent(out) :: hi, lo
      real(dp) :: t

      t = SPLITTER * a
      hi = t - (t - a)
      lo = a - hi
   end subroutine split

   pure subroutine two_prod(a, b, p, e)
      !! p = fl(a*b) and e = the exact discarded error, so p + e == a * b.
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: p, e
      real(dp) :: ahi, alo, bhi, blo

      p = a * b
      call split(a, ahi, alo)
      call split(b, bhi, blo)
      e = ((ahi * bhi - p) + ahi * blo + alo * bhi) + alo * blo
   end subroutine two_prod

   pure type(dd) function dd_add(a, b) result(r)
      !! Sum of two double-doubles, renormalised.
      type(dd), intent(in) :: a, b
      real(dp) :: s1, s2, t1, t2, u1, u2

      ! Every `two_sum` writes to variables it does not also read. Passing the
      ! same name as both an intent(in) and an intent(out) argument is aliasing,
      ! which Fortran leaves undefined -- and the first version of this routine
      ! did exactly that, renormalising through `two_sum(s1, s2, s1, s2)`. It
      ! compiled without a murmur and returned a sum of roughly zero.
      call two_sum(a%hi, b%hi, s1, s2)
      call two_sum(a%lo, b%lo, t1, t2)
      s2 = s2 + t1
      call two_sum(s1, s2, u1, u2)
      u2 = u2 + t2
      call two_sum(u1, u2, r%hi, r%lo)
   end function dd_add

   pure type(dd) function dd_mul(a, b) result(r)
      !! Product of two double-doubles, renormalised.
      type(dd), intent(in) :: a, b
      real(dp) :: p1, p2

      call two_prod(a%hi, b%hi, p1, p2)
      p2 = p2 + a%hi * b%lo + a%lo * b%hi
      call two_sum(p1, p2, r%hi, r%lo)
   end function dd_mul

   pure type(dd) function dd_neg(a) result(r)
      type(dd), intent(in) :: a
      r%hi = -a%hi; r%lo = -a%lo
   end function dd_neg

   pure type(dd) function dd_sub(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = dd_add(a, dd_neg(b))
   end function dd_sub

   pure type(dd) function dd_div(a, b) result(r)
      !! Long division: a double-precision quotient, then one correction pass
      !! against the exact remainder.
      !!
      !! The remainder is what makes it work. `a - q1*b` is computed in
      !! double-double, so the leading terms cancel exactly rather than to
      !! within a rounding, and what survives is the part the first quotient
      !! missed.
      type(dd), intent(in) :: a, b
      real(dp) :: q1, q2
      type(dd) :: r1, prod

      q1 = a%hi / b%hi
      prod = dd_mul(b, dd_from(q1))
      r1 = dd_sub(a, prod)
      q2 = (r1%hi + r1%lo) / b%hi
      call two_sum(q1, q2, r%hi, r%lo)
   end function dd_div

   pure type(dd) function dd_sqrt(a) result(r)
      !! One Newton step on a double-precision square root.
      !!
      !! x1 = x0 + (a - x0^2) / (2 x0) doubles the correct digits, so a
      !! 53-bit seed lands at the ~106 the type carries. The residual `a - x0^2`
      !! is formed in double-double for the same reason as in the division:
      !! computed in double it would be pure rounding noise.
      type(dd), intent(in) :: a
      real(dp) :: x0
      type(dd) :: resid, corr

      if (a%hi <= 0.0_dp) then
         r%hi = 0.0_dp; r%lo = 0.0_dp
         return
      end if
      x0 = sqrt(a%hi)
      resid = dd_sub(a, dd_mul(dd_from(x0), dd_from(x0)))
      corr = dd_div(resid, dd_from(2.0_dp * x0))
      r = dd_add(dd_from(x0), corr)
   end function dd_sqrt

   pure type(dd) function dd_abs(a) result(r)
      type(dd), intent(in) :: a
      if (a%hi < 0.0_dp) then
         r = dd_neg(a)
      else
         r = a
      end if
   end function dd_abs

   ! Mixed dd/real64 forms, so a literal in the shared bodies needs no wrapper.
   pure type(dd) function dd_add_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_add(a, dd_from(b))
   end function dd_add_r

   pure type(dd) function dd_r_add(a, b) result(r)
      real(dp), intent(in) :: a
      type(dd), intent(in) :: b
      r = dd_add(dd_from(a), b)
   end function dd_r_add

   pure type(dd) function dd_sub_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_sub(a, dd_from(b))
   end function dd_sub_r

   pure type(dd) function dd_r_sub(a, b) result(r)
      real(dp), intent(in) :: a
      type(dd), intent(in) :: b
      r = dd_sub(dd_from(a), b)
   end function dd_r_sub

   pure type(dd) function dd_mul_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_mul(a, dd_from(b))
   end function dd_mul_r

   pure type(dd) function dd_r_mul(a, b) result(r)
      real(dp), intent(in) :: a
      type(dd), intent(in) :: b
      r = dd_mul(dd_from(a), b)
   end function dd_r_mul

   pure type(dd) function dd_div_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_div(a, dd_from(b))
   end function dd_div_r

   pure type(dd) function dd_r_div(a, b) result(r)
      real(dp), intent(in) :: a
      type(dd), intent(in) :: b
      r = dd_div(dd_from(a), b)
   end function dd_r_div

   ! Comparisons go through the full value, not just `hi`: two numbers can share
   ! a leading double and differ in the tail, which is the entire reason for the
   ! type.
   pure logical function dd_lt(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = (a%hi < b%hi) .or. (a%hi == b%hi .and. a%lo < b%lo)
   end function dd_lt

   pure logical function dd_lt_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_lt(a, dd_from(b))
   end function dd_lt_r

   pure logical function dd_gt(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = (a%hi > b%hi) .or. (a%hi == b%hi .and. a%lo > b%lo)
   end function dd_gt

   pure logical function dd_gt_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_gt(a, dd_from(b))
   end function dd_gt_r

   pure logical function dd_le(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = .not. dd_gt(a, b)
   end function dd_le

   pure logical function dd_le_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_le(a, dd_from(b))
   end function dd_le_r

   pure logical function dd_ge(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = .not. dd_lt(a, b)
   end function dd_ge

   pure logical function dd_ge_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_ge(a, dd_from(b))
   end function dd_ge_r

   pure logical function dd_eq(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = (a%hi == b%hi) .and. (a%lo == b%lo)
   end function dd_eq

   pure logical function dd_eq_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_eq(a, dd_from(b))
   end function dd_eq_r

   pure type(dd) function dd_exp(a) result(r)
      !! exp for double-double, by range reduction and a Taylor series.
      !!
      !! Two reductions. First `a = k ln2 + t` puts `t` in [-ln2/2, ln2/2]; then
      !! `t` is halved NSQ more times and the result squared back up.
      !!
      !! **NSQ is a trade, and it was measured rather than guessed.** Each
      !! squaring doubles the relative error, so NSQ of 16 puts a floor of
      !! 2^16 * 1e-32 under the answer -- 2.3e-27, measured. Fewer squarings
      !! mean a larger series argument and more terms. Over [-35, 35]: NSQ 16
      !! gives 2.3e-27, NSQ 8 gives 1.4e-29, NSQ 6 gives 3.3e-30. Below 6 the
      !! twelve terms below stop converging and it gets worse again.
      type(dd), intent(in) :: a
      integer, parameter :: NSQ = 6
      type(dd) :: t, term, sum, sq
      real(dp) :: k
      integer :: i

      if (a%hi < -700.0_dp) then
         r = dd_from(0.0_dp)
         return
      end if

      k = anint(dd_to_dp(dd_mul(a, DD_INV_LN2)))
      t = dd_sub(a, dd_mul(DD_LN2, dd_from(k)))
      ! Halve NSQ times: exp(t) = exp(t/2^NSQ)^(2^NSQ)
      t = dd_mul(t, dd_from(2.0_dp**(-NSQ)))

      ! 1 + t + t^2/2! + ...
      sum = dd_add(dd_from(1.0_dp), t)
      term = t
      do i = 2, 12
         term = dd_mul(term, t)
         term = dd_div(term, dd_from(real(i, dp)))
         sum = dd_add(sum, term)
      end do

      do i = 1, NSQ
         sum = dd_mul(sum, sum)
      end do

      r = dd_mul(sum, dd_from(2.0_dp**int(k)))
   end function dd_exp

   pure type(dd) function dd_erf(a) result(r)
      !! erf for double-double.
      !!
      !! Below the crossover this uses the *confluent* series
      !!
      !!     erf(x) = (2x/sqrt(pi)) exp(-x^2) sum_n (2x^2)^n / (1.3.5...(2n+1))
      !!
      !! rather than the textbook alternating one. Every term is positive, so
      !! there is no cancellation; the alternating series loses roughly x^2/ln10
      !! digits to it and is useless here by x = 4, which is inside the range
      !! the Boys function asks for.
      type(dd), intent(in) :: a
      type(dd) :: x, x2, term, sum, two_x2
      real(dp) :: den
      integer :: n

      x = dd_abs(a)
      if (x%hi > 6.5_dp) then
         r = dd_from(1.0_dp)
         if (a%hi < 0.0_dp) r = dd_neg(r)
         return
      end if

      x2 = dd_mul(x, x)
      two_x2 = dd_mul(x2, dd_from(2.0_dp))
      term = dd_from(1.0_dp)
      sum = dd_from(1.0_dp)
      do n = 1, 200
         den = real(2*n + 1, dp)
         term = dd_div(dd_mul(term, two_x2), dd_from(den))
         sum = dd_add(sum, term)
         if (abs(term%hi) < 1.0e-36_dp * abs(sum%hi)) exit
      end do

      r = dd_mul(dd_mul(DD_2_SQRTPI, x), dd_mul(dd_exp(dd_neg(x2)), sum))
      if (a%hi < 0.0_dp) r = dd_neg(r)
   end function dd_erf

   ! dd_erfc is deliberately absent: Jorge has a bitwise-reproducible erfc that
   ! should be the reference here, and guessing at an asymptotic series before
   ! seeing it would mean writing something that has to agree with libcint bit
   ! for bit without ever having been compared to it.
   !
   ! The two call sites that need it are both `erfc(a) - erfc(b)` -- in
   ! cint_fmt_body.inc and cint_wheeler_body.inc -- so whatever lands has to be
   ! accurate in the tail, where `1 - erf(x)` is not.

end module cint_dd
