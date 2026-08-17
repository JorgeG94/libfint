!
! D8 acceptance: the libcint_fortran module, as a caller sees it.
!
! Only built with WITH_FORTRAN_BACKEND, where libcint_fortran is the port and
! libcint_interface is still the C -- so one program can put the two side by
! side through exactly the API metalquicha uses, 2D (SLOTS, *) tables and all.
!
! This is a different test from int1e_check/int2e_check/int3c2e_check, which
! call the port's own entry points directly.  What it adds is the shim: the
! sequence association that reinterprets the caller's tables, the env extent
! it works out for itself, the output dimensions it derives rather than being
! told.  Any of those could be wrong while every engine test still passes.
!
! THE BAR IS BIT-IDENTICAL UP TO FIVE RYS ROOTS AND 1e-11 ABOVE, and the split
! is the one D3 predicted rather than one chosen to make this pass.  Up to five
! roots the quadrature comes from the closed forms, which the port reproduces
! exactly; above that it comes from a tridiagonal eigensolve, and D3 replaced
! the C's vendored copy of DSTEMR with LAPACK's.  Asked to referee the two by
! a symmetry the integral itself obeys -- (ij|kl) = (ji|kl), which the exact
! answer satisfies and neither engine quite does -- the port comes out very
! slightly ahead:
!
!     roots   quartets       C          port
!       4        3720     1.13e-14    1.13e-14      (identical)
!       5        2528     3.17e-14    3.17e-14      (identical)
!       6         714     8.25e-14    7.43e-14
!       7          64     2.02e-13    1.92e-13
!
! so this is the accuracy floor of two different eigensolvers, not a defect on
! either side.  The tolerance sits two orders above that floor.
!
program backend_check
   use iso_c_binding
   use libcint_interface, only: cint1e_ovlp_cart, cint1e_ovlp_sph, &
                                cint1e_kin_cart,  cint1e_kin_sph,  &
                                cint1e_nuc_cart,  cint1e_nuc_sph,  &
                                cint2e_cart,      cint2e_sph,      &
                                cint3c2e_cart,    cint3c2e_sph,    &
                                cint2c2e_cart,    cint2c2e_sph,    &
                                cint1e_ipnuc_cart,  cint1e_ipnuc_sph,  &
                                cint1e_iprinv_cart, cint1e_iprinv_sph, &
                                cint2e_ip1_cart,    cint2e_ip1_sph,    &
                                cint3c2e_ip1_cart,  cint3c2e_ip1_sph,  &
                                cint3c2e_ip2_cart,  cint3c2e_ip2_sph,  &
                                cint2c2e_ip1_cart,  cint2c2e_ip1_sph,  &
                                CINTcgto_cart,    CINTcgto_spheric, &
                                ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, ANG_OF
   use libcint_fortran
   use cint_test_systems
   implicit none

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, l, shls(4), ret
   real(c_double), allocatable :: cbuf(:), fbuf(:)
   integer :: ibasis, sph, di, dj, dk, dl, nq, stride, nc, nf
   integer :: ncmp, nbad, nexact, nloose
   real(dp) :: worst
   character(len=80) :: wat
   type(c_ptr) :: opt

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   allocate(cbuf(400000), fbuf(400000))
   ncmp = 0; nbad = 0; nexact = 0; nloose = 0; worst = 0; wat = "(none)"

   do ibasis = 1, n_reference_basis
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ibasis, bas, env, nbas, off)

      ! The bookkeeping first: if these disagree the buffers are the wrong
      ! size and every comparison below is meaningless.
      do i = 0, nbas-1
         call check_int(int(CINTcgto_cart(i, bas)), int(libcint_cgto_cart(i, bas)), &
                        "cgto_cart", ibasis, i)
         call check_int(int(CINTcgto_spheric(i, bas)), int(libcint_cgto_sph(i, bas)), &
                        "cgto_sph", ibasis, i)
      end do

      ! An optimizer, requested and then ignored, because that is what a
      ! caller does and the null it gets back has to be harmless.
      call libcint_2e_sph_optimizer(opt, atm, 8_ip, bas, nbas, env)

      stride = max(1, int(nbas)/6)

      do sph = 0, 1
         ! ---- one electron, every pair ----
         do i = 0, nbas-1
         do j = 0, nbas-1
            call dims2(sph, i, j, bas, di, dj)
            nq = di*dj
            shls(1)=i; shls(2)=j
            if (sph == 1) then
               ret = cint1e_ovlp_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ovlp_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ovlp_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ovlp_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, nq, "ovlp", ibasis, sph, i, j, 0, 0)

            if (sph == 1) then
               ret = cint1e_kin_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_kin_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_kin_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_kin_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, nq, "kin", ibasis, sph, i, j, 0, 0)

            if (sph == 1) then
               ret = cint1e_nuc_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_nuc_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_nuc_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_nuc_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, nq, "nuc", ibasis, sph, i, j, 0, 0)

            ! ---- two centre, every pair ----
            call dims2(sph, i, j, bas, di, dj)
            if (sph == 1) then
               ret = cint2c2e_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_2c2e_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint2c2e_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_2c2e_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj, "2c2e", ibasis, sph, i, j, 0, 0)

            ! ---- the derivative families, which were stubs here until the
            ! port grew them.  Three components each.
            call dims2(sph, i, j, bas, di, dj)
            if (sph == 1) then
               ret = cint1e_ipnuc_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipnuc_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_ipnuc_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_ipnuc_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*3, "ipnuc", ibasis, sph, i, j, 0, 0, 1)

            if (sph == 1) then
               ret = cint1e_iprinv_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_iprinv_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint1e_iprinv_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               ret = libcint_1e_iprinv_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*3, "iprinv", ibasis, sph, i, j, 0, 0, 1)

            if (sph == 1) then
               ret = cint2c2e_ip1_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_2c2e_ip1_sph(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            else
               ret = cint2c2e_ip1_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_2c2e_ip1_cart(fbuf, shls(1:2), atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*3, "2c2e_ip1", ibasis, sph, i, j, 0, 0, 1)
         end do
         end do

         ! ---- three and four centre, strided ----
         do i = 0, nbas-1, stride
         do j = 0, nbas-1, stride
         do k = 0, nbas-1, stride
            call dims2(sph, i, j, bas, di, dj)
            call dims1(sph, k, bas, dk)
            shls(1)=i; shls(2)=j; shls(3)=k; shls(4)=0
            if (sph == 1) then
               ret = cint3c2e_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
            else
               ret = cint3c2e_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*dk, "3c2e", ibasis, sph, i, j, k, 0)

            if (sph == 1) then
               ret = cint3c2e_ip1_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_ip1_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
            else
               ret = cint3c2e_ip1_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_ip1_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*dk*3, "3c2e_ip1", ibasis, sph, i, j, k, 0, 1)

            if (sph == 1) then
               ret = cint3c2e_ip2_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_ip2_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
            else
               ret = cint3c2e_ip2_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
               ret = libcint_3c2e_ip2_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
            end if
            call cmp(cbuf, fbuf, di*dj*dk*3, "3c2e_ip2", ibasis, sph, i, j, k, 0, 1)

            do l = 0, nbas-1, stride
               call dims1(sph, l, bas, dl)
               shls(4)=l
               if (sph == 1) then
                  ret = cint2e_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_sph(fbuf, shls, atm, 8_ip, bas, nbas, env, opt)
               else
                  ret = cint2e_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_cart(fbuf, shls, atm, 8_ip, bas, nbas, env, opt)
               end if
               call cmp(cbuf, fbuf, di*dj*dk*dl, "2e", ibasis, sph, i, j, k, l)

               if (sph == 1) then
                  ret = cint2e_ip1_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ip1_sph(fbuf, shls, atm, 8_ip, bas, nbas, env)
               else
                  ret = cint2e_ip1_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
                  ret = libcint_2e_ip1_cart(fbuf, shls, atm, 8_ip, bas, nbas, env)
               end if
               call cmp(cbuf, fbuf, di*dj*dk*dl*3, "2e_ip1", ibasis, sph, i, j, k, l, 1)
            end do
         end do
         end do
         end do
      end do

      call libcint_del_optimizer(opt)
   end do

   print '(A,I0)',         "  values compared      : ", ncmp
   print '(A,ES10.2)',     "  worst diff / block   : ", worst
   print '(A,I0)',         "  blocks above 5 roots : ", nloose
   print '(A,I0)',         "  inexact at <=5 roots : ", nexact
   print '(A,I0)',         "  over tolerance       : ", nbad
   if (nbad > 0 .or. nexact > 0) then
      print '(A,A)',       "  worst at             : ", trim(wat)
      stop 1
   end if
   print '(A)', "  RESULT: PASS through libcint_fortran"

contains

   subroutine dims1(sp, sh, b, d)
      integer, intent(in) :: sp
      integer(c_int), intent(in) :: sh, b(BAS_SLOTS, *)
      integer, intent(out) :: d
      if (sp == 1) then
         d = int(CINTcgto_spheric(sh, b))
      else
         d = int(CINTcgto_cart(sh, b))
      end if
   end subroutine dims1

   subroutine dims2(sp, sa, sb, b, da, db)
      integer, intent(in) :: sp
      integer(c_int), intent(in) :: sa, sb, b(BAS_SLOTS, *)
      integer, intent(out) :: da, db
      call dims1(sp, sa, b, da)
      call dims1(sp, sb, b, db)
   end subroutine dims2

   subroutine check_int(a, b, what, ib, sh)
      integer, intent(in) :: a, b, ib, sh
      character(len=*), intent(in) :: what
      ncmp = ncmp + 1
      if (a /= b) then
         nbad = nbad + 1
         if (worst < 1.0_dp) then
            worst = 1.0_dp
            write(wat,'(A,A,I0,A,I0,A,I0,A,I0)') trim(what), " basis=", ib, &
               " shell=", sh, " C=", a, " F=", b
         end if
      end if
   end subroutine check_int

   ! Rys order for a set of shells: the expression the envs setup uses.
   !
   ! `inc` is the integral's own ng increment total.  It is not optional
   ! bookkeeping: a derivative raises the angular ceiling on the shell it acts
   ! on, so int2e_ip1 over four p shells needs six roots where int2e over the
   ! same four needs five.  Counting from the bare l undercounts, and the
   ! block then gets judged at tol=0 when it went through the eigensolver --
   ! which is where the port and the C are allowed to differ.  That mistake
   ! is what 169,353 "failures" at a worst of 2.4e-14 looked like.
   integer function rys_order(b, s1, s2, s3, s4, inc)
      integer(c_int), intent(in) :: b(BAS_SLOTS, *)
      integer, intent(in) :: s1, s2, s3, s4
      integer, intent(in), optional :: inc
      integer :: extra
      extra = 0
      if (present(inc)) extra = inc
      rys_order = (int(b(ANG_OF, s1+1)) + int(b(ANG_OF, s2+1)) &
                 + int(b(ANG_OF, s3+1)) + int(b(ANG_OF, s4+1)) + extra) / 2 + 1
   end function rys_order

   subroutine cmp(c, f, n, what, ib, sp, ii, jj, kk, ll, inc)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, sp, ii, jj, kk, ll
      character(len=*), intent(in) :: what
      integer, intent(in), optional :: inc
      integer :: m, nr
      real(dp) :: r, scal, tol
      if (n <= 0) return
      nr = rys_order(bas, ii, jj, kk, ll, inc)
      if (nr <= 5) then
         tol = 0.0_dp
      else
         tol = 1.0e-11_dp
         nloose = nloose + 1
      end if
      ! Scale by the largest integral in the block rather than element by
      ! element: an entry near zero divided by itself measures nothing.
      scal = maxval(abs(c(1:n)))
      if (scal <= 0.0_dp) scal = 1.0_dp
      do m = 1, n
         r = abs(f(m)-c(m)) / scal
         ncmp = ncmp + 1
         if (r > worst) then
            worst = r
            write(wat,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(what), &
               " basis=", ib, " sph=", sp, " shls (",ii,",",jj,",",kk,",",ll, &
               ") m=", m, " roots=", nr
         end if
         if (r > tol) nbad = nbad + 1
         if (r /= 0.0_dp .and. nr <= 5) nexact = nexact + 1
      end do
   end subroutine cmp

end program backend_check
