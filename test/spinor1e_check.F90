!
! D9: the one-electron spinor integrals, end to end, against the C.
!
! Sweeps every shell pair of all five reference basis sets.  The reference
! systems leave KAPPA_OF at zero, which is the general case -- both j = l+1/2
! and j = l-1/2 present, 4l+2 components -- so the sweep is also rerun with
! kappa forced to each sign in turn, because those select the other
! coefficient block and a wrong choice there is invisible at kappa = 0.
!
! int1e_spsp is here for the D9 acceptance, and because it is the generator's
! first 1e spinor entry point: same description, same gout, same primitive
! loop as _cart and _sph, and only the transform on the way out differs.
!
! The bar is bit-identical.  int1e_ovlp_spinor has no Rys quadrature in it at
! all; int1e_nuc_spinor does, but at (l+l)/2+1 roots it stays inside the
! closed forms for every basis here, and so does int1e_spsp at (l+1+l+1)/2+1.
!
program spinor1e_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, &
                                KAPPA_OF, ANG_OF, CINTcgto_spinor
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_1e_spinor, only: int1e_ovlp_spinor, int1e_nuc_spinor
   use cint_gen_intor3, only: int1e_spsp_spinor
   implicit none

   interface
      function c_ovlp(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int1e_ovlp_spinor')
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_ovlp
      end function
      function c_spsp(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int1e_spsp_spinor')
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_spsp
      end function
      function c_nuc(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int1e_nuc_spinor')
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_nuc
      end function
   end interface

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, shls(2), ret
   complex(c_double_complex), allocatable :: cbuf(:), fbuf(:)
   integer :: ibasis, kmode, di, dj, nq, dims(0:3), fshls(0:1)
   integer :: ncmp, nbad, l
   character(len=64) :: wat
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ncmp = 0; nbad = 0; wat = "(none)"

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
         l = int(bas(ANG_OF, i))
         select case (kmode)
         case (1); bas(KAPPA_OF, i) = int(-l - 1, c_int)
         case (2); if (l > 0) bas(KAPPA_OF, i) = int(l, c_int)
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

      do i = 0, nbas-1
      do j = 0, nbas-1
         di = int(CINTcgto_spinor(i, bas)); dj = int(CINTcgto_spinor(j, bas))
         if (di <= 0 .or. dj <= 0) cycle
         nq = di*dj
         allocate(cbuf(nq), fbuf(nq))
         cbuf = 0; fbuf = 0
         shls(1) = i; shls(2) = j
         fshls(0) = i; fshls(1) = j
         dims(0) = di; dims(1) = dj; dims(2) = 1; dims(3) = 1

         ret = c_ovlp(cbuf, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                      c_null_ptr, c_null_ptr)
         hv = int1e_ovlp_spinor(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cbuf, fbuf, nq, "ovlp", ibasis, kmode, i, j)

         cbuf = 0; fbuf = 0
         ret = c_nuc(cbuf, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                     c_null_ptr, c_null_ptr)
         hv = int1e_nuc_spinor(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cbuf, fbuf, nq, "nuc", ibasis, kmode, i, j)

         cbuf = 0; fbuf = 0
         ret = c_spsp(cbuf, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                      c_null_ptr, c_null_ptr)
         hv = int1e_spsp_spinor(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmp(cbuf, fbuf, nq, "spsp", ibasis, kmode, i, j)

         deallocate(cbuf, fbuf)
      end do
      end do
      deallocate(fatm, fbas, fenv)
   end do
   end do

   print '(A,I0)', "  values compared : ", ncmp
   print '(A,I0)', "  differing       : ", nbad
   if (nbad > 0) then
      print '(A,A)', "  first at        : ", trim(wat)
      stop 1
   end if
   print '(A)', "  RESULT: PASS (bit-identical)"

contains

   subroutine cmp(c, f, n, what, ib, km, ii, jj)
      complex(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, km, ii, jj
      character(len=*), intent(in) :: what
      integer :: m
      do m = 1, n
         ncmp = ncmp + 1
         if (c(m) /= f(m)) then
            nbad = nbad + 1
            if (nbad == 1) write(wat,'(A,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(what), &
               " basis=", ib, " kmode=", km, " shls (", ii, ",", jj, ") m=", m
         end if
      end do
   end subroutine cmp

end program spinor1e_check
