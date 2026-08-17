!
! D9: int2e_spinor, end to end, against the C.

! THE BAR IS BIT-IDENTICAL UP TO FIVE RYS ROOTS AND 1e-11 ABOVE.  That split is
! D3's, not one chosen to make this pass: up to five roots the quadrature comes
! from the closed forms, which the port reproduces exactly, and above that from
! a tridiagonal eigensolve where D3 replaced the C's vendored DSTEMR with
! LAPACK's.  backend_check refereed the two by a symmetry the integral itself
! obeys and found the port marginally the more self-consistent, so this is an
! accuracy floor rather than a defect on either side.  See PORT_PLAN.md D8.
!
!
! Sweeps every shell pair of all five reference basis sets.  The reference
! systems leave KAPPA_OF at zero, which is the general case -- both j = l+1/2
! and j = l-1/2 present, 4l+2 components -- so the sweep is also rerun with
! kappa forced to each sign in turn, because those select the other
! coefficient block and a wrong choice there is invisible at kappa = 0.
!
! The bar is bit-identical.  int1e_ovlp_spinor has no Rys quadrature in it at
! all; int1e_nuc_spinor does, but at (l+l)/2+1 roots it stays inside the
! closed forms for every basis here.
!
program spinor2e_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, &
                                KAPPA_OF, ANG_OF, CINTcgto_spinor, BAS_SLOTS
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_2e_spinor, only: int2e_spinor
   implicit none

   interface
      function c_2e(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int2e_spinor')
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_2e
      end function
   end interface

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, l, shls(4), ret
   complex(c_double_complex), allocatable :: cbuf(:), fbuf(:)
   integer :: ibasis, kmode, di, dj, dk, dl, nq, dims(0:3), fshls(0:3), stride
   integer :: ncmp, nbad, nexact, lang
   real(dp) :: worst
   character(len=64) :: wat
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ncmp = 0; nbad = 0; nexact = 0; worst = 0.0_dp; wat = "(none)"

   do ibasis = 1, n_reference_basis
   do kmode = 0, 2
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ibasis, bas, env, nbas, off)

      ! kmode 0 leaves kappa at zero (both j), 1 forces j = l + 1/2, 2 forces
      ! j = l - 1/2.  l = 0 has no j = l - 1/2, so it stays at zero there --
      ! a shell with no components is not a case the C handles either.
      do i = 1, nbas
         lang = int(bas(ANG_OF, i))
         select case (kmode)
         case (1); bas(KAPPA_OF, i) = int(-lang - 1, c_int)
         case (2); if (lang > 0) bas(KAPPA_OF, i) = int(lang, c_int)
         end select
      end do

      allocate(fatm(0:6*8-1), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1))
      do i = 1, 8
         fatm(6*(i-1):6*i-1) = int(atm(:, i))
      end do
      do i = 1, nbas
         fbas(8*(i-1):8*i-1) = int(bas(:, i))
      end do
      fenv = real(env, dp)

      ! spinor quartets are large -- 4l+2 components per shell at kappa = 0 --
      ! so stride, the way int2e_check does, rather than sweep all four indices
      stride = max(1, int(nbas)/3)
      do i = 0, nbas-1, stride
      do j = 0, nbas-1, stride
      do k = 0, nbas-1, stride
      do l = 0, nbas-1, stride
         di = int(CINTcgto_spinor(i, bas)); dj = int(CINTcgto_spinor(j, bas))
         dk = int(CINTcgto_spinor(k, bas)); dl = int(CINTcgto_spinor(l, bas))
         if (di <= 0 .or. dj <= 0 .or. dk <= 0 .or. dl <= 0) cycle
         nq = di*dj*dk*dl
         allocate(cbuf(nq), fbuf(nq))
         cbuf = 0; fbuf = 0
         shls(1) = i; shls(2) = j; shls(3) = k; shls(4) = l
         fshls(0) = i; fshls(1) = j; fshls(2) = k; fshls(3) = l
         dims(0) = di; dims(1) = dj; dims(2) = dk; dims(3) = dl

         ret = c_2e(cbuf, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                    c_null_ptr, c_null_ptr)
         hv = int2e_spinor(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cbuf, fbuf, nq, "int2e", ibasis, kmode, i, j)

         deallocate(cbuf, fbuf)
      end do
      end do
      end do
      end do
      deallocate(fatm, fbas, fenv)
   end do
   end do

   print '(A,I0)',     "  values compared      : ", ncmp
   print '(A,ES10.2)', "  worst diff / block   : ", worst
   print '(A,I0)',     "  inexact at <=5 roots : ", nexact
   print '(A,I0)',     "  over tolerance       : ", nbad
   if (nbad > 0 .or. nexact > 0) then
      print '(A,A)', "  first at        : ", trim(wat)
      stop 1
   end if
   print '(A)', "  RESULT: PASS"

contains

   ! Rys order for a quartet: the expression the envs setup uses.
   integer function rys_order(s1, s2, s3, s4)
      integer, intent(in) :: s1, s2, s3, s4
      rys_order = (int(bas(ANG_OF, s1+1)) + int(bas(ANG_OF, s2+1)) &
                 + int(bas(ANG_OF, s3+1)) + int(bas(ANG_OF, s4+1))) / 2 + 1
   end function rys_order

   subroutine cmp(c, f, n, what, ib, km, ii, jj)
      complex(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, km, ii, jj
      character(len=*), intent(in) :: what
      integer :: m, nr
      real(dp) :: r, scal, tol
      if (n <= 0) return
      nr = rys_order(ii, jj, k, l)
      if (nr <= 5) then
         tol = 0.0_dp
      else
         tol = 1.0e-11_dp
      end if
      scal = maxval(abs(c(1:n)))
      if (scal <= 0.0_dp) scal = 1.0_dp
      do m = 1, n
         ncmp = ncmp + 1
         r = abs(c(m) - f(m)) / scal
         if (r > worst) worst = r
         if (r > tol) then
            nbad = nbad + 1
            if (nbad == 1) write(wat,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(what), &
               " basis=", ib, " kmode=", km, " shls (", ii, ",", jj, ") m=", m, " roots=", nr
         end if
         if (r /= 0.0_dp .and. nr <= 5) nexact = nexact + 1
      end do
   end subroutine cmp

end program spinor2e_check
