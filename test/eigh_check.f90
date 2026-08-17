! The whole eigensolver against eigh.c, compared BITWISE.
!
! Bitwise rather than by difference, for two reasons.  It is the actual claim
! -- the port reproduces the C's doubles, not merely close ones.  And the
! reference itself returns NaN for a few wildly scaled inputs, where |a-b| is
! NaN however equal the two are; comparing the representations says "identical"
! for those, which is the truth.
program eigh_check
   use iso_fortran_env, only: int64
   use cint_const, only: dp
   use cint_eigh_dstemr, only: dstemr_diagonalize
   implicit none
   integer, parameter :: MX = 32
   real(dp) :: d(0:MX-1), e(0:MX-1), w(0:MX-1), z(0:MX*MX-1)
   real(dp) :: cw(0:MX-1), cz(0:MX*MX-1)
   integer  :: n, i, info, cinfo, ncase, nbad, nnan, u, ios, ndiff
   character(len=256) :: path
   if (command_argument_count() < 1) then
      print *, "usage: eigh_check <reference.bin>"
      stop 2
   end if
   call get_command_argument(1, path)
   open(newunit=u, file=path, access="stream", form="unformatted", status="old")
   ncase = 0; nbad = 0; nnan = 0
   do
      read(u, iostat=ios) n
      if (ios /= 0) exit
      read(u) d(0:n-1); read(u) e(0:n-1)
      read(u) cinfo;    read(u) cw(0:n-1); read(u) cz(0:n*n-1)
      w = 0.0_dp; z = 0.0_dp
      info = dstemr_diagonalize(n, d, e, w, z)
      ncase = ncase + 1
      ndiff = 0
      do i = 0, n-1
         if (differs(w(i), cw(i))) ndiff = ndiff + 1
      end do
      do i = 0, n*n-1
         if (differs(z(i), cz(i))) ndiff = ndiff + 1
      end do
      if (any(isnan_(cw(0:n-1)))) nnan = nnan + 1
      if (ndiff /= 0 .or. info /= cinfo) then
         nbad = nbad + 1
         if (nbad <= 4) print '(A,I3,A,I0,A,I0,A,I0)', "  n=", n, "  differing words: ", &
            ndiff, "   info ", info, " vs ", cinfo
      end if
   end do
   close(u)
   print '(A,I0,A,I0,A,I0)', "  cases: ", ncase, "   (reference NaN in ", nnan, ")"
   print '(A,I0)',           "  differing: ", nbad
   if (nbad > 0) stop 1
   print *, " OK: bit-identical to eigh.c"
contains
   ! Two NaNs count as the same answer.  A NaN payload is not specified by
   ! IEEE-754 and neither implementation chooses one deliberately; what
   ! matters is that both said "no answer" for the same input.
   logical function differs(a, b)
      real(dp), intent(in) :: a, b
      if (isnan_(a) .and. isnan_(b)) then
         differs = .false.
      else
         differs = bits(a) /= bits(b)
      end if
   end function differs

   integer(int64) function bits(x)
      real(dp), intent(in) :: x
      bits = transfer(x, 1_int64)
   end function bits
   elemental logical function isnan_(x)
      real(dp), intent(in) :: x
      isnan_ = (x /= x)
   end function isnan_
end program eigh_check
