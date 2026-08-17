!
! D8: the density-fitting integrals, (ij|k) and (i|k), against the C.
!
! Same shape as int2e_check: stride through the shell triples (and pairs) of
! all five reference basis sets, Cartesian and spherical, and require
! bit-identity.
!
! There is no separate auxiliary basis here.  The C's own test suite does the
! same thing -- it hands int3c2e the ordinary orbital shells for all three
! indices -- and it is the harder case, because the k shell then reaches the
! same angular momenta as i and j rather than the s and p an auxiliary set is
! mostly made of.
!
program int3c2e_check
   use iso_c_binding
   use libcint_interface
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_3c2e, only: int3c2e_cart, int3c2e_sph, int2c2e_cart, int2c2e_sph
   implicit none
   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, di, dj, dk, shls(4), ret
   real(c_double), allocatable :: cbuf(:), fbuf(:)
   integer :: ibasis, sph, dims(0:3), ncmp3, nbad3, ncmp2, nbad2, nq, stride
   real(dp) :: worst3, worst2
   character(len=72) :: wat3, wat2
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:), fshls(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ncmp3 = 0; nbad3 = 0; worst3 = 0; wat3 = "(none)"
   ncmp2 = 0; nbad2 = 0; worst2 = 0; wat2 = "(none)"

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

      ! three indices instead of four, so a finer stride than int2e_check's
      ! costs the same and covers more angular-momentum combinations
      stride = max(1, int(nbas)/15)

      ! ---- three centre ----
      do sph = 0, 1
      do i = 0, nbas-1, stride
      do j = 0, nbas-1, stride
      do k = 0, nbas-1, stride
         if (sph == 1) then
            di = CINTcgto_spheric(i,bas); dj = CINTcgto_spheric(j,bas)
            dk = CINTcgto_spheric(k,bas)
         else
            di = CINTcgto_cart(i,bas); dj = CINTcgto_cart(j,bas)
            dk = CINTcgto_cart(k,bas)
         end if
         nq = di*dj*dk
         allocate(cbuf(nq), fbuf(nq))
         cbuf = 0; fbuf = 0
         shls(1)=i; shls(2)=j; shls(3)=k
         fshls(0)=i; fshls(1)=j; fshls(2)=k
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         if (sph == 1) then
            ret = cint3c2e_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
            hv = int3c2e_sph(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         else
            ret = cint3c2e_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
            hv = int3c2e_cart(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         end if
         call cmp(cbuf, fbuf, nq, ibasis, sph, i, j, k, ncmp3, nbad3, worst3, wat3)
         deallocate(cbuf, fbuf)
      end do
      end do
      end do
      end do

      ! ---- two centre ----
      do sph = 0, 1
      do i = 0, nbas-1
      do k = 0, nbas-1
         if (sph == 1) then
            di = CINTcgto_spheric(i,bas); dk = CINTcgto_spheric(k,bas)
         else
            di = CINTcgto_cart(i,bas); dk = CINTcgto_cart(k,bas)
         end if
         nq = di*dk
         allocate(cbuf(nq), fbuf(nq))
         cbuf = 0; fbuf = 0
         shls(1)=i; shls(2)=k
         fshls(0)=i; fshls(1)=k
         dims(0)=di; dims(1)=dk; dims(2)=1; dims(3)=1
         if (sph == 1) then
            ret = cint2c2e_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
            hv = int2c2e_sph(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         else
            ret = cint2c2e_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env, c_null_ptr)
            hv = int2c2e_cart(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         end if
         call cmp(cbuf, fbuf, nq, ibasis, sph, i, 0, k, ncmp2, nbad2, worst2, wat2)
         deallocate(cbuf, fbuf)
      end do
      end do
      end do

      deallocate(fatm, fbas, fenv, fshls)
   end do

   print '(A,I0,A,ES10.2,A,I0)', "  int3c2e  values=", ncmp3, &
         "  worst=", worst3, "  differing=", nbad3
   if (nbad3 > 0) print '(A,A)', "           at ", trim(wat3)
   print '(A,I0,A,ES10.2,A,I0)', "  int2c2e  values=", ncmp2, &
         "  worst=", worst2, "  differing=", nbad2
   if (nbad2 > 0) print '(A,A)', "           at ", trim(wat2)
   if (nbad3 > 0 .or. nbad2 > 0) stop 1
   print '(A)', "  RESULT: PASS (bit-identical)"
contains
   subroutine cmp(c, f, n, ib, sp, ii, jj, kk, ncmp, nbad, worst, wat)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, sp, ii, jj, kk
      integer, intent(inout) :: ncmp, nbad
      real(dp), intent(inout) :: worst
      character(len=*), intent(inout) :: wat
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
            write(wat,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') "basis=",ib," sph=",sp, &
               " shls (",ii,",",jj,",",kk,") m=",m
         end if
         if (r /= 0.0_dp) nbad = nbad + 1
      end do
   end subroutine
end program int3c2e_check
