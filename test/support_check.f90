!
! Acceptance test for the support layer (D4).
!
!   support_check <c2s-dump> <blas-dump>
!
! Everything here is expected to be BIT-IDENTICAL to the C.  Unlike the Rys
! layer, nothing in D4 changes an algorithm: the Cartesian-to-spherical
! transforms are the same sparse sums in the same order, gto_norm is the same
! expression, and a transpose is a permutation.  So the bar is exact equality
! and any difference at all is a failure.
!
program support_check
   use cint_const,    only: dp
   use cint_bas,      only: cint_len_cart, cint_gto_norm
   use cint_cart2sph, only: cint_c2s_bra_sph, cint_c2s_ket_sph, RESULT_IN_GCART
   use cint_blas,     only: cint_dmat_transpose, cint_dplus_transpose, cint_daxpy2v
   implicit none

   character(len=256) :: p1, p2
   integer  :: nbad1, ncmp1, nbad2, ncmp2
   real(dp) :: worst1
   character(len=64) :: wat1

   if (command_argument_count() < 2) then
      print '(A)', "usage: support_check <c2s-dump> <blas-dump>"
      stop 2
   end if
   call get_command_argument(1, p1)
   call get_command_argument(2, p2)

   call run_c2s(trim(p1), nbad1, ncmp1, worst1, wat1)
   print '(A)',            "  [1] cart-to-spherical transforms and gto_norm"
   print '(A,I0)',         "      values compared : ", ncmp1
   print '(A,ES10.2,A,A)', "      worst rel diff  : ", worst1, "   at ", trim(wat1)
   print '(A,I0)',         "      differing       : ", nbad1

   call run_blas(trim(p2), nbad2, ncmp2)
   print '(A)',            "  [2] transposes and axpy"
   print '(A,I0)',         "      values compared : ", ncmp2
   print '(A,I0)',         "      differing       : ", nbad2

   print '(A)', ""
   if (nbad1 == 0 .and. nbad2 == 0) then
      print '(A)', "  RESULT: PASS (bit-identical)"
   else
      print '(A)', "  RESULT: FAIL"
      stop 1
   end if

contains

   subroutine run_c2s(path, nbad, ncmp, worst, wat)
      character(*), intent(in)  :: path
      integer,      intent(out) :: nbad, ncmp
      real(dp),     intent(out) :: worst
      character(*), intent(out) :: wat

      integer  :: u, ios, kind_, l, n, nf, nd, lds, i, j, loc
      real(dp) :: gcart(0:4095), cres(0:4095), gsph(0:4095), a, cv, fv, rel

      open(newunit=u, file=path, access='stream', form='unformatted', status='old')
      nbad = 0; ncmp = 0; worst = 0.0_dp; wat = "(none)"
      do
         read(u, iostat=ios) kind_
         if (ios /= 0) exit
         read(u) l, n
         nf = cint_len_cart(l)
         nd = 2*l + 1

         if (kind_ == 2) then
            read(u) a, cv
            fv = cint_gto_norm(l, a)
            call note(reldiff(fv, cv), "gto_norm l=", l, n, nbad, ncmp, worst, wat)
            cycle
         end if

         read(u) gcart(0:nf*n-1)
         if (kind_ == 0) then
            read(u) cres(0:nd*n-1)
            gsph = 0.0_dp
            loc = cint_c2s_bra_sph(gsph, n, gcart, l)
            do i = 0, nd*n - 1
               if (loc == RESULT_IN_GCART) then
                  rel = reldiff(gcart(i), cres(i))
               else
                  rel = reldiff(gsph(i), cres(i))
               end if
               call note(rel, "bra l=", l, n, nbad, ncmp, worst, wat)
            end do
         else
            lds = n + 2
            read(u) cres(0:nd*lds-1)
            gsph = 0.0_dp
            loc = cint_c2s_ket_sph(gsph, gcart, lds, n, l)
            do i = 0, nd - 1
               do j = 0, n - 1
                  if (loc == RESULT_IN_GCART) then
                     rel = reldiff(gcart(i*n+j), cres(i*n+j))
                  else
                     rel = reldiff(gsph(i*lds+j), cres(i*lds+j))
                  end if
                  call note(rel, "ket l=", l, n, nbad, ncmp, worst, wat)
               end do
            end do
         end if
      end do
      close(u)
   end subroutine run_c2s

   subroutine note(r, tag, l, n, nbad, ncmp, worst, wat)
      real(dp),     intent(in)    :: r
      character(*), intent(in)    :: tag
      integer,      intent(in)    :: l, n
      integer,      intent(inout) :: nbad, ncmp
      real(dp),     intent(inout) :: worst
      character(*), intent(inout) :: wat

      ncmp = ncmp + 1
      if (r > worst) then
         worst = r
         write(wat,'(A,I0,A,I0)') tag, l, " n=", n
      end if
      if (r /= 0.0_dp) nbad = nbad + 1
   end subroutine note

   subroutine run_blas(path, nbad, ncmp)
      character(*), intent(in)  :: path
      integer,      intent(out) :: nbad, ncmp

      integer  :: u, ios, kind_, m, n, i
      real(dp) :: a(0:4095), ct(0:4095), ft(0:4095)
      real(dp) :: y(0:4095), cv(0:4095), fv(0:4095)

      open(newunit=u, file=path, access='stream', form='unformatted', status='old')
      nbad = 0; ncmp = 0
      do
         read(u, iostat=ios) kind_
         if (ios /= 0) exit
         read(u) m, n

         if (kind_ == 2) then
            read(u) a(0:m-1), y(0:m-1), cv(0:m-1)
            call cint_daxpy2v(m, 1.37_dp, a, y, fv)
            do i = 0, m - 1
               ncmp = ncmp + 1
               if (fv(i) /= cv(i)) nbad = nbad + 1
            end do
            cycle
         end if

         read(u) a(0:m*n-1), ct(0:m*n-1)
         ! same seed the C used, so dplus_transpose accumulates onto the
         ! same starting values
         do i = 0, m*n - 1
            ft(i) = (mod(i,5) + 1) * 0.11_dp
         end do
         if (kind_ == 0) then
            call cint_dmat_transpose(ft, a, m, n)
         else
            call cint_dplus_transpose(ft, a, m, n)
         end if
         do i = 0, m*n - 1
            ncmp = ncmp + 1
            if (ft(i) /= ct(i)) nbad = nbad + 1
         end do
      end do
      close(u)
   end subroutine run_blas

   pure real(dp) function reldiff(a, b)
      real(dp), intent(in) :: a, b
      if (b /= 0.0_dp) then
         reldiff = abs(a - b) / abs(b)
      else
         reldiff = abs(a - b)
      end if
   end function reldiff

end program support_check
