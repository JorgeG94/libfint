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
! **This breaks under -ffast-math.**  Every one of these routines depends on the
! compiler *not* reassociating `(a + b) - a`, which is the whole mechanism.
! -Ofast, -ffast-math and nvfortran's default -fast will silently turn these
! into zero error terms and leave double precision wearing a double-double type.
! The build must keep IEEE semantics for this file.
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

end module cint_dd
