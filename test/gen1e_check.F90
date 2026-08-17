!
! D6 acceptance: the GENERATED one-electron kernels against the C.
!
! Every shell pair of the five reference basis sets, Cartesian and spherical,
! for four integrals chosen to cover the distinct shapes the emitter has to
! handle: no operation (ovlp), a derivative on the ket (kin), an
! origin-shifted operator (r), and a derivative on the bra (ipovlp).
!
! The bar is bit-identical.  D1 established that "it compiles" proves nothing
! here -- the first version of that spike emitted code which compiled and was
! wrong -- so this compares values, and getting it to pass turned up three
! real defects: a missing dispatch through envs%f_gout, a missing origin
! vector for the RC operations, and an associativity difference between C's
! `+=` and restating the target in Fortran.
!
program gen1e_check
   use iso_c_binding
   use libcint_interface
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_gen_intor1
   use cint_1e, only: int1e_ovlp_cart, int1e_ovlp_sph
   use cint_gen_grad1, only: int1e_ipovlp_cart, int1e_ipovlp_sph
   implicit none
   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, di, dj, shls(2), ret
   real(c_double), allocatable :: cbuf(:), fbuf(:)
   integer :: ibasis, sph, dims(0:1), which, ncmp(4), nbad(4), nc
   real(dp) :: rel, worst(4)
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:), fshls(:)
   real(dp), allocatable :: fenv(:)
   character(len=14) :: nm(4)
   logical :: hv

   ! int1e_r is not in libcint_interface, so bind it here rather than widen
   ! that module for a test.
   interface
      function cint1e_r_cart(buf, shls, atm, natm, bas, nbas, env) &
            bind(C, name='cint1e_r_cart') result(r)
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in)  :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         real(c_double), intent(in)  :: env(*)
         integer(c_int) :: r
      end function
      function cint1e_r_sph(buf, shls, atm, natm, bas, nbas, env) &
            bind(C, name='cint1e_r_sph') result(r)
         import :: c_double, c_int
         real(c_double), intent(out) :: buf(*)
         integer(c_int), intent(in)  :: shls(*), atm(*), bas(*)
         integer(c_int), value, intent(in) :: natm, nbas
         real(c_double), intent(in)  :: env(*)
         integer(c_int) :: r
      end function
   end interface
   nm = [character(len=14)::"int1e_ovlp","int1e_kin","int1e_r","int1e_ipovlp"]

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ncmp = 0; nbad = 0; worst = 0

   do ibasis = 1, n_reference_basis
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ibasis, bas, env, nbas, off)
      allocate(fatm(0:6*8-1), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1), fshls(0:1))
      do i = 1, 8
         fatm(6*(i-1):6*i-1) = int(atm(:, i))
      end do
      do i = 1, nbas
         fbas(8*(i-1):8*i-1) = int(bas(:, i))
      end do
      fenv = real(env, dp)

      do which = 1, 4
      do sph = 0, 1
      do i = 0, nbas-1
      do j = 0, nbas-1
         if (sph == 1) then
            di = CINTcgto_spheric(i, bas); dj = CINTcgto_spheric(j, bas)
         else
            di = CINTcgto_cart(i, bas);    dj = CINTcgto_cart(j, bas)
         end if
         nc = di*dj
         if (which == 3) nc = nc*3       ! int1e_r has 3 components
         if (which == 4) nc = nc*3       ! int1e_ipovlp has 3
         allocate(cbuf(nc), fbuf(nc))
         cbuf = 0; fbuf = 0
         shls(1) = i; shls(2) = j
         fshls(0) = i; fshls(1) = j
         dims(0) = di; dims(1) = dj
         select case (which)
         case (1)
            if (sph==1) then
               ret = cint1e_ovlp_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               hv = int1e_ovlp_sph(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            else
               ret = cint1e_ovlp_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               hv = int1e_ovlp_cart(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            end if
         case (2)
            if (sph==1) then
               ret = cint1e_kin_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               hv = int1e_kin_sph(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            else
               ret = cint1e_kin_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               hv = int1e_kin_cart(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            end if
         case (3)
            if (sph==1) then
               ret = cint1e_r_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               hv = int1e_r_sph(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            else
               ret = cint1e_r_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               hv = int1e_r_cart(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            end if
         case (4)
            if (sph==1) then
               ret = cint1e_ipovlp_sph(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               hv = int1e_ipovlp_sph(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            else
               ret = cint1e_ipovlp_cart(cbuf, shls, atm, 8_c_int, bas, nbas, env)
               hv = int1e_ipovlp_cart(fbuf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            end if
         end select
         call cmp(cbuf, fbuf, nc, which)
         deallocate(cbuf, fbuf)
      end do
      end do
      end do
      end do
      deallocate(fatm, fbas, fenv, fshls)
   end do

   do which = 1, 4
      print '(A,A,A,I0,A,ES10.2,A,I0)', "  ", nm(which), " values=", ncmp(which), &
         "  worst=", worst(which), "  differing=", nbad(which)
   end do
   if (sum(nbad) > 0) stop 1
   print '(A)', "  RESULT: PASS (bit-identical)"
contains
   subroutine cmp(c, f, n, w)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, w
      integer :: k
      real(dp) :: r
      do k = 1, n
         if (c(k) /= 0.0_dp) then
            r = abs(f(k)-c(k))/abs(c(k))
         else
            r = abs(f(k)-c(k))
         end if
         ncmp(w) = ncmp(w) + 1
         worst(w) = max(worst(w), r)
         if (r /= 0.0_dp) nbad(w) = nbad(w) + 1
      end do
   end subroutine
end program gen1e_check
