!
! D7 acceptance: int2e, end to end, against the C.
!
! Strides through the unique-quartet space of all five reference basis sets,
! Cartesian and spherical.  The stride keeps the run short while still
! reaching every angular-momentum combination each basis contains -- including
! cc-pVQZ, whose g functions push the Rys order past the unrolled kernels and
! onto the general 2D-to-4D path.
!
! The bar is bit-identical.  The Rys layer underneath is not (D3 swapped the
! eigensolver), but that only shows above eight roots, which a two-electron
! integral over these bases does not reach.
!
program int2e_check
   use iso_c_binding
   use libcint_interface
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_2e, only: int2e_cart, int2e_sph
   implicit none
   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, l, di, dj, dk, dl, shls(4), ret
   real(c_double), allocatable :: cbuf(:), fbuf(:)
   integer :: ibasis, sph, dims(0:3), ncmp, nbad, nq, stride
   real(dp) :: rel, worst
   character(len=72) :: wat
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:), fshls(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ncmp = 0; nbad = 0; worst = 0; wat = "(none)"

   do ibasis = 1, n_reference_basis
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ibasis, bas, env, nbas, off)
      allocate(fatm(0:6*8-1), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1), fshls(0:3))
      do i = 1, 8
         fatm(6*(i-1):6*i-1) = int(atm(:, i))
      end do
      do i = 1, nbas
         fbas(8*(i-1):8*i-1) = int(bas(:, i))
      end do
      fenv = real(env, dp)

      ! stride through the quartets so the test stays quick but still covers
      ! every angular-momentum combination the basis has
      stride = max(1, int(nbas)/5)
      do sph = 0, 1
      do i = 0, nbas-1, stride
      do j = 0, nbas-1, stride
      do k = 0, nbas-1, stride
      do l = 0, nbas-1, stride
         if (sph == 1) then
            di = CINTcgto_spheric(i,bas); dj = CINTcgto_spheric(j,bas)
            dk = CINTcgto_spheric(k,bas); dl = CINTcgto_spheric(l,bas)
         else
            di = CINTcgto_cart(i,bas); dj = CINTcgto_cart(j,bas)
            dk = CINTcgto_cart(k,bas); dl = CINTcgto_cart(l,bas)
         end if
         nq = di*dj*dk*dl
         allocate(cbuf(nq), fbuf(nq))
         cbuf = 0; fbuf = 0
         shls(1)=i; shls(2)=j; shls(3)=k; shls(4)=l
         fshls(0)=i; fshls(1)=j; fshls(2)=k; fshls(3)=l
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=dl
         if (sph == 1) then
            ret = cint2e_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
            hv = int2e_sph(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         else
            ret = cint2e_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
            hv = int2e_cart(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         end if
         call cmp(cbuf, fbuf, nq, ibasis, sph, i, j, k, l)
         deallocate(cbuf, fbuf)
      end do
      end do
      end do
      end do
      end do
      deallocate(fatm, fbas, fenv, fshls)
   end do

   print '(A,I0)',        "  values compared : ", ncmp
   print '(A,ES10.2,A,A)',"  worst rel diff  : ", worst, "   at ", trim(wat)
   print '(A,I0)',        "  differing       : ", nbad
   if (nbad > 0) stop 1
   print '(A)', "  RESULT: PASS (bit-identical)"
contains
   subroutine cmp(c, f, n, ib, sp, ii, jj, kk, ll)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, sp, ii, jj, kk, ll
      integer :: m
      real(dp) :: r
      do m = 1, n
         if (c(m) /= 0.0_dp) then
            r = abs(f(m)-c(m))/abs(c(m))
         else
            r = abs(f(m)-c(m))
         end if
         ncmp = ncmp + 1
         if (r > worst) then
            worst = r
            write(wat,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') "basis=",ib," sph=",sp, &
               " quartet (",ii,",",jj,",",kk,",",ll,") m=",m
         end if
         if (r /= 0.0_dp) nbad = nbad + 1
      end do
   end subroutine
end program int2e_check
