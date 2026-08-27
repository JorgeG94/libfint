program q_check
   !! The shim must be bit-identical to real128 where IEEE-754 says there is
   !! only one right answer, and within a measured ulp budget where it does not.
   !!
   !! The standard mandates correct rounding for + - * / and sqrt, and for
   !! nothing else. Those five therefore have a single correct result that
   !! every implementation produces, so anything short of bit-equality there is
   !! lossy plumbing -- which is what this test was written to catch, and still
   !! catches.
   !!
   !! The transcendentals carry no such guarantee, and the compilers duly
   !! disagree. gfortran's real128 intrinsics resolve to glibc's *f128 family,
   !! the same symbols the shim calls, so they agree to the bit; ifx resolves
   !! them into Intel's libimf -- __expq, __erfq, __erfcq -- which is a
   !! different implementation of a function nobody is obliged to round
   !! correctly. Over 20,000 points that is 7,442 erfc values differing, 126
   !! erf and 1 exp, never by more than 2 ulp, while sqrt -- mandated -- agrees
   !! on every one. Demanding bit-equality on erfc would be demanding the shim
   !! reproduce whichever libm the compiler happened to link, which no shim can
   !! do everywhere.
   !!
   !! Compared as raw bits, all 128 of them, on both sides of the split:
   !! converting to double first would forgive 60 low bits of drift, and doing
   !! so has already produced a false agreement twice in this work.
   use iso_fortran_env, only: real64, real128, int64
   use cint_quad, only: quad, quad_from, operator(+), operator(-), &
                        operator(*), operator(/), sqrt, exp, erf, erfc
   implicit none
   integer, parameter :: dp = real64, qp = real128
   integer, parameter :: npt = 2000
   !! Worst measured gap between the shim and the local real128 is 2 ulp, on
   !! ifx's erfc, and 0 on gfortran, where both sides are the same glibc.
   !! Doubling that leaves room for another libm without leaving room for a
   !! plumbing bug, which would show up as a gap of many bits, not four.
   integer(int64), parameter :: ulp_budget = 4_int64
   integer :: i, k
   integer :: nbad(6)
   integer(int64) :: worst(6)
   real(dp) :: x
   real(qp) :: ref
   type(quad) :: z
   logical :: failed
   character(len=12) :: nm(6)
   !! .true. where IEEE-754 mandates correct rounding, and the budget is zero.
   logical, parameter :: mandated(6) = [.true., .true., .true., .false., .false., .false.]
   nm = [character(len=12) :: "sqrt", "x*x - 1", "(x+1)/3", "exp", "erf", "erfc"]

   nbad = 0
   worst = 0

   do i = 1, npt
      x = real(i, dp)*0.0031_dp
      z = sqrt(quad_from(x)); ref = sqrt(real(x, qp)); call acc(1, z, ref)
      z = quad_from(x)*quad_from(x) - quad_from(1.0_dp)
      ref = real(x, qp)*real(x, qp) - 1.0_qp
      call acc(2, z, ref)
      z = (quad_from(x) + quad_from(1.0_dp))/quad_from(3.0_dp)
      ref = (real(x, qp) + 1.0_qp)/3.0_qp
      call acc(3, z, ref)
      z = exp(quad_from(x)); ref = exp(real(x, qp)); call acc(4, z, ref)
      z = erf(quad_from(x)); ref = erf(real(x, qp)); call acc(5, z, ref)
      z = erfc(quad_from(x)); ref = erfc(real(x, qp)); call acc(6, z, ref)
   end do

   failed = .false.
   do k = 1, 6
      print '(a,a12,a,i6,a,i0,a,i0,a)', "  ", nm(k), "  differing: ", nbad(k), &
         " of ", npt, "   worst gap: ", worst(k), " ulp"
      if (mandated(k)) then
         if (nbad(k) /= 0) failed = .true.
      else
         if (worst(k) > ulp_budget) failed = .true.
      end if
   end do

   if (failed) error stop "C shim does not track real128"
   print *, "OK -- exact on the IEEE-mandated operations, within budget elsewhere"

contains

   subroutine acc(k, a, b)
      integer, intent(in) :: k
      type(quad), intent(in) :: a
      real(qp), intent(in) :: b
      integer(int64) :: d
      if (mandated(k)) then
         !! Raw bits, not magnitudes: an exact operation must reproduce the
         !! sign too, and x*x - 1 is negative over most of this sweep.
         if (any(transfer(a, [0_int64, 0_int64]) /= transfer(b, [0_int64, 0_int64]))) then
            nbad(k) = nbad(k) + 1
            worst(k) = max(worst(k), ulp_gap(a, b))
         end if
      else
         d = ulp_gap(a, b)
         if (d /= 0_int64) nbad(k) = nbad(k) + 1
         worst(k) = max(worst(k), d)
      end if
   end subroutine acc

   integer(int64) function ulp_gap(a, b) result(g)
      !! Distance in units in the last place, taken on the 128-bit pattern.
      !!
      !! For finite values of one sign the IEEE encoding is monotone read as an
      !! integer, so the difference of the two patterns IS the number of
      !! representable values between them -- across exponent boundaries
      !! included, which is what makes this the right measure for a function
      !! whose result spans 300 decades. Sign is masked off and the subtraction
      !! done unsigned over both words; anything that will not fit in an int64
      !! is reported as huge, since a gap that size is a bug and not a rounding
      !! difference.
      type(quad), intent(in) :: a
      real(qp), intent(in) :: b
      integer(int64) :: wa(2), wb(2), hi, lo, t
      wa = transfer(a, [0_int64, 0_int64])
      wb = transfer(b, [0_int64, 0_int64])
      wa(2) = iand(wa(2), int(z'7FFFFFFFFFFFFFFF', int64))
      wb(2) = iand(wb(2), int(z'7FFFFFFFFFFFFFFF', int64))
      if (merge(ult(wa(1), wb(1)), ult(wa(2), wb(2)), wa(2) == wb(2))) then
         t = wa(1); wa(1) = wb(1); wb(1) = t
         t = wa(2); wa(2) = wb(2); wb(2) = t
      end if
      lo = wa(1) - wb(1)
      hi = wa(2) - wb(2)
      if (ult(wa(1), wb(1))) hi = hi - 1_int64
      if (hi /= 0_int64 .or. lo < 0_int64) then
         g = huge(0_int64)
      else
         g = lo
      end if
   end function ulp_gap

   logical function ult(x1, y1) result(v)
      !! Unsigned <, on a pattern Fortran only has a signed integer for.
      integer(int64), intent(in) :: x1, y1
      integer(int64), parameter :: flip = int(z'8000000000000000', int64)
      v = ieor(x1, flip) < ieor(y1, flip)
   end function ult

end program q_check
