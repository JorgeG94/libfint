!
! Checks the Rys root finder and the Schmidt path against the C (D3).
!
!   rys_core_check <roots-dump> <schmidt-dump>
!
! Both Fortran precisions are expected to be BIT-IDENTICAL to their own C
! counterpart -- the double path to CINTrys_schmidt, the quadruple path to
! CINTqrys_schmidt, since both of those are binary128.  The C's long double is
! reported alongside for information only: it is 80-bit x87 and differs for
! that reason, not because either side is wrong.
!
! Matching bit for bit required reproducing two places where the C's own three
! copies have drifted apart; see cint_schmidt.F90.
!
program rys_core_check
   use cint_const, only: dp
   use cint_quad, only: quad, quad_from, quad_to_dp, operator(+), operator(-), &
                      operator(*), operator(/), operator(<), operator(>), &
                      operator(<=), operator(>=), operator(==), operator(**), &
                      assignment(=), sqrt, abs, exp, erf, erfc, max, min
   use cint_find_roots, only: cint_polynomial_roots
   use cint_schmidt_dp, only: schmidt_dp => rys_schmidt
   use cint_schmidt_qp, only: schmidt_qp => rys_schmidt
   implicit none

   character(len=256) :: p1, p2
   integer :: nfail

   call get_command_argument(1, p1)
   call get_command_argument(2, p2)
   nfail = 0
   call check_roots(trim(p1), nfail)
   call check_schmidt(trim(p2), nfail)

   print '(A)', ""
   if (nfail == 0) then
      print '(A)', "  RESULT: PASS"
   else
      print '(A,I0,A)', "  RESULT: FAIL (", nfail, ")"
      stop 1
   end if

contains

   subroutine check_roots(path, nfail)
      character(*), intent(in) :: path
      integer, intent(inout) :: nfail
      integer :: u, ios, nroots, scenario, cerr, ferr, nroots1, i, n, nbad
      real(dp) :: cs(0:32*32-1), croots(0:31), froots(0:31)
      real(dp) :: rel, worst

      open(newunit=u, file=path, access='stream', form='unformatted', status='old')
      n = 0; nbad = 0; worst = 0.0_dp
      do
         read(u, iostat=ios) nroots
         if (ios /= 0) exit
         read(u) scenario
         read(u) cerr
         nroots1 = nroots + 1
         read(u) cs(0:nroots1*nroots1-1)
         read(u) croots(0:nroots-1)
         froots = 0.0_dp
         ferr = cint_polynomial_roots(froots, cs, nroots)
         if (ferr /= cerr) then
            nbad = nbad + 1
            cycle
         end if
         if (cerr /= 0) cycle
         do i = 0, nroots-1
            rel = rd(froots(i), croots(i))
            worst = max(worst, rel)
            if (rel /= 0.0_dp) nbad = nbad + 1
            n = n + 1
         end do
      end do
      close(u)
      print '(A)',        "  [1] polynomial roots, dp against the C"
      print '(A,I0)',     "      roots compared : ", n
      print '(A,ES10.2)', "      worst rel diff : ", worst
      if (nbad > 0) then
         print '(A,I0)',  "      NOT bit-identical, differing: ", nbad
         nfail = nfail + 1
      end if
   end subroutine check_roots

   subroutine check_schmidt(path, nfail)
      character(*), intent(in) :: path
      integer, intent(inout) :: nfail
      integer :: u, ios, n, variant, cerr, ferr, i
      real(dp) :: x, lo, cu(0:63), cw(0:63), fu(0:63), fw(0:63)
      real(dp) :: rel, worst(0:2)
      integer  :: ncmp(0:2), nbad(0:2), nerr(0:2)
      character(len=13) :: nm(0:2)
      nm = [character(len=13) :: "dp  vs C dp  ", "qp  vs C ld  ", "qp  vs C quad"]

      open(newunit=u, file=path, access='stream', form='unformatted', status='old')
      worst = 0.0_dp; ncmp = 0; nbad = 0; nerr = 0
      do
         read(u, iostat=ios) n
         if (ios /= 0) exit
         read(u) variant, cerr, x, lo, cu(0:n-1), cw(0:n-1)
         fu = 0.0_dp; fw = 0.0_dp
         if (variant == 0) then
            ferr = schmidt_dp(n, x, lo, fu, fw)
         else
            ferr = schmidt_qp(n, x, lo, fu, fw)
         end if
         if (ferr /= cerr) then
            nerr(variant) = nerr(variant) + 1
            cycle
         end if
         if (cerr /= 0) cycle
         do i = 0, n-1
            rel = max(rd(fu(i), cu(i)), rd(fw(i), cw(i)))
            worst(variant) = max(worst(variant), rel)
            if (rel /= 0.0_dp) nbad(variant) = nbad(variant) + 1
            ncmp(variant) = ncmp(variant) + 1
         end do
      end do
      close(u)

      print '(A)', "  [2] Schmidt roots and weights"
      do i = 0, 2
         print '(A,A,A,I0,A,ES9.2,A,I0,A,I0)', "      ", nm(i), " n=", ncmp(i), &
            "  worst=", worst(i), "  differing=", nbad(i), "  err-mismatch=", nerr(i)
      end do
      ! A count of zero for the quad variant is not a pass by accident: the
      ! reference omits it when libcint was built without quadmath, which is
      ! what happens under a compiler that does not ship it -- icx, for one.
      ! Say so, rather than printing a bare zero that reads like a bug.
      if (ncmp(2) == 0) then
         print '(A)', "      quad: not compared -- the reference libcint has" // &
                      " no quadmath, so there is"
         print '(A)', "            no binary128 C path to compare against."
      end if
      ! Variant 1 is the C's 80-bit path; it is not expected to agree.
      if (nbad(0) > 0 .or. nerr(0) > 0) nfail = nfail + 1
      if (nbad(2) > 0 .or. nerr(2) > 0) nfail = nfail + 1
   end subroutine check_schmidt

   real(dp) function rd(a, b)
      real(dp), intent(in) :: a, b
      if (b /= 0.0_dp) then
         rd = abs(a-b)/abs(b)
      else
         rd = abs(a-b)
      end if
   end function rd

end program rys_core_check
