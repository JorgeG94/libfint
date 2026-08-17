!
! D5 first light: int1e_ovlp, end to end, from the Fortran port.
!
! Runs every shell pair of the five C2H6 reference basis sets through both
! implementations, Cartesian and spherical, and compares.  This is the
! milestone test for the whole 1e architecture -- it exercises the envs
! derived type, the procedure-pointer gout dispatch, the workspace that
! replaced MALLOC_INSTACK, primitive-pair screening, the Obara-Saika
! recursion, the contraction step, and both cart2sph paths.
!
! The bar is bit-identical.  Nothing in this layer changes an algorithm or an
! order of operations, so anything else would be a mistake -- and unlike the
! Rys layer there is no eigensolver substitution to excuse a difference.
!
! Note the two index conventions in play, which is what D0 caught a bug in.
! The libcint_interface module the C side uses has 1-based slot constants
! because it indexes a Fortran array column; the port keeps the C's flat,
! 0-based tables.  The conversion is explicit below.
!
program int1e_check
   use iso_c_binding
   use libcint_interface
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_1e, only: int1e_ovlp_cart, int1e_ovlp_sph
   implicit none
   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, di, dj, shls(2), ret
   real(c_double), allocatable :: cbuf(:), fbuf(:)
   integer :: ibasis, sph, dims(0:1), ncmp, nbad
   real(dp) :: rel, worst
   character(len=64) :: wat
   type(cint_ws) :: ws
   ! flat 0-based copies for the Fortran port, which keeps the C's layout
   integer, allocatable :: fatm(:), fbas(:), fshls(:)
   real(dp), allocatable :: fenv(:)

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ncmp = 0; nbad = 0; worst = 0; wat = "(none)"

   do ibasis = 1, n_reference_basis
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ibasis, bas, env, nbas, off)

      ! the port takes flat, 0-based tables with 0-based slot constants
      allocate(fatm(0:6*8-1), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1), fshls(0:1))
      do i = 1, 8
         fatm(6*(i-1):6*i-1) = int(atm(:, i))
      end do
      do i = 1, nbas
         fbas(8*(i-1):8*i-1) = int(bas(:, i))
      end do
      fenv = real(env, dp)

      do sph = 0, 1
      do i = 0, nbas-1
         do j = 0, nbas-1
            if (sph == 1) then
               di = CINTcgto_spheric(i, bas); dj = CINTcgto_spheric(j, bas)
            else
               di = CINTcgto_cart(i, bas);    dj = CINTcgto_cart(j, bas)
            end if
            allocate(cbuf(di*dj), fbuf(di*dj))
            cbuf = 0; fbuf = 0
            shls(1) = i; shls(2) = j
            if (sph == 1) then
               ret = cint1e_ovlp_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
            else
               ret = cint1e_ovlp_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
            end if
            fshls(0) = i; fshls(1) = j
            dims(0) = di; dims(1) = dj
            if (sph == 1) then
               ret = merge(1, 0, int1e_ovlp_sph(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws))
            else
               ret = merge(1, 0, int1e_ovlp_cart(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws))
            end if
            call cmp(cbuf, fbuf, di*dj, ibasis, sph, i, j)
            deallocate(cbuf, fbuf)
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
   subroutine cmp(c, f, n, ib, sp, ii, jj)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, sp, ii, jj
      integer :: k
      real(dp) :: r
      do k = 1, n
         if (c(k) /= 0.0_dp) then
            r = abs(f(k)-c(k))/abs(c(k))
         else
            r = abs(f(k)-c(k))
         end if
         ncmp = ncmp + 1
         if (r > worst) then
            worst = r
            write(wat,'(A,I0,A,I0,A,I0,A,I0,A,I0)') "basis=",ib," sph=",sp, &
               " shells (",ii,",",jj,") k=",k
         end if
         if (r /= 0.0_dp) nbad = nbad + 1
      end do
   end subroutine
end program int1e_check
