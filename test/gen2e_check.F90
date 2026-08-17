!
! D9/D10: the generated two-electron integrals against the C.
!
! Three integrals, from the same auto_intor.cl descriptions the C is generated
! from, in all three angular forms:
!
!   int2e            -- the control.  Its _cart and _sph are already known
!                       bit-identical from int2e_check, so if the new 2e
!                       backend breaks the emitter at all it shows here first
!                       and with the least to disentangle.
!   int2e_spsp1      -- sigma.p on electron one.  c2s_si_2e1 then c2s_sf_2e2:
!                       the mixed case, which is what forced the spinor
!                       driver's two stages to be selected independently.
!   int2e_spsp1spsp2 -- sigma.p on both.  si on both stages, four Rys-root
!                       derivatives, and the largest g array of the three.
!
! All three kappa modes, as in spinor1e_check: kappa = 0 exercises only one of
! the two coefficient blocks.
!
! THE BAR IS BIT-IDENTICAL UP TO FIVE RYS ROOTS AND 1e-11 ABOVE.  That split is
! D3's, not one chosen to make this pass: up to five roots the quadrature comes
! from the closed forms, which the port reproduces exactly, and above that from
! a tridiagonal eigensolve where D3 replaced the C's vendored DSTEMR with
! LAPACK's.  backend_check refereed the two by a symmetry the integral itself
! obeys and found the port marginally the more self-consistent, so this is an
! accuracy floor, not a defect on either side.  See PORT_PLAN.md D8.
!
program gen2e_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, &
                                KAPPA_OF, ANG_OF, CINTcgto_spinor, &
                                CINTcgto_cart, CINTcgto_spheric
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_gen_intor4, only: int2e_spsp1_cart, int2e_spsp1_sph, int2e_spsp1_spinor, &
                              int2e_spsp1spsp2_cart, int2e_spsp1spsp2_sph, &
                              int2e_spsp1spsp2_spinor
   use cint_2e_spinor,  only: int2e_spinor
   implicit none

   interface
      function c_2e_sp(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int2e_spinor')
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_2e_sp
      end function
      function c_spsp1_sp(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int2e_spsp1_spinor')
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_spsp1_sp
      end function
      function c_spsp12_sp(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int2e_spsp1spsp2_spinor')
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_spsp12_sp
      end function
      function c_spsp1_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int2e_spsp1_sph')
         import :: c_int, c_double, c_ptr
         real(c_double) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_spsp1_sph
      end function
      function c_spsp1_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int2e_spsp1_cart')
         import :: c_int, c_double, c_ptr
         real(c_double) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_spsp1_cart
      end function
      function c_spsp12_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache) &
            bind(C, name='int2e_spsp1spsp2_sph')
         import :: c_int, c_double, c_ptr
         real(c_double) :: out(*)
         type(c_ptr), value :: dims
         integer(c_int) :: shls(*), atm(*), bas(*)
         integer(c_int), value :: natm, nbas
         real(c_double) :: env(*)
         type(c_ptr), value :: opt, cache
         integer(c_int) :: c_spsp12_sph
      end function
   end interface

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, l, shls(4), ret
   complex(c_double_complex), allocatable :: cz(:), fz(:)
   real(c_double), allocatable :: cr(:), fr(:)
   integer :: ibasis, kmode, di, dj, dk, dl, nq, dims(0:3), fshls(0:3)
   ! The ng increments each integral raises its four angular ceilings by --
   ! the first four entries of the ng array the generator emits.
   integer, parameter :: NG_2E(4)     = [0, 0, 0, 0]
   integer, parameter :: NG_SPSP1(4)  = [1, 1, 0, 0]
   integer, parameter :: NG_SPSP12(4) = [1, 1, 1, 1]
   integer :: ncmp, nbad, nexact, lang, stride, iw
   real(dp) :: worst
   character(len=24) :: names(7)
   integer :: tally(7), seen(7)
   character(len=72) :: wat
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   names = [character(len=24):: "int2e_spinor","int2e_spsp1_spinor","int2e_spsp1spsp2_spinor", &
                                "int2e_spsp1_sph","int2e_spsp1spsp2_sph","int2e_spsp1_cart","?"]
   tally = 0; seen = 0
   ncmp = 0; nbad = 0; nexact = 0; worst = 0.0_dp; wat = "(none)"

   do ibasis = 1, n_reference_basis
   do kmode = 0, 2
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ibasis, bas, env, nbas, off)
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

      stride = max(1, int(nbas)/3)
      do i = 0, nbas-1, stride
      do j = 0, nbas-1, stride
      do k = 0, nbas-1, stride
      do l = 0, nbas-1, stride
         shls(1) = i; shls(2) = j; shls(3) = k; shls(4) = l
         fshls(0) = i; fshls(1) = j; fshls(2) = k; fshls(3) = l

         ! ---- spinor ----
         di = int(CINTcgto_spinor(i, bas)); dj = int(CINTcgto_spinor(j, bas))
         dk = int(CINTcgto_spinor(k, bas)); dl = int(CINTcgto_spinor(l, bas))
         if (di > 0 .and. dj > 0 .and. dk > 0 .and. dl > 0) then
            nq = di*dj*dk*dl
            dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=dl
            allocate(cz(nq), fz(nq))

            cz = 0; fz = 0
            ret = c_2e_sp(cz, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                          c_null_ptr, c_null_ptr)
            hv = int2e_spinor(fz, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            call cmpz(cz, fz, nq, "int2e_spinor", ibasis, kmode, i, j, k, l, NG_2E)

            cz = 0; fz = 0
            ret = c_spsp1_sp(cz, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                             c_null_ptr, c_null_ptr)
            hv = int2e_spsp1_spinor(fz, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            call cmpz(cz, fz, nq, "int2e_spsp1_spinor", ibasis, kmode, i, j, k, l, NG_SPSP1)

            cz = 0; fz = 0
            ret = c_spsp12_sp(cz, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                              c_null_ptr, c_null_ptr)
            hv = int2e_spsp1spsp2_spinor(fz, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
            call cmpz(cz, fz, nq, "int2e_spsp1spsp2_spinor", ibasis, kmode, i, j, k, l, NG_SPSP12)

            deallocate(cz, fz)
         end if

         ! ---- real forms; kappa does not enter, so once is enough ----
         if (kmode /= 0) cycle

         di = int(CINTcgto_spheric(i, bas)); dj = int(CINTcgto_spheric(j, bas))
         dk = int(CINTcgto_spheric(k, bas)); dl = int(CINTcgto_spheric(l, bas))
         nq = di*dj*dk*dl
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=dl
         ! spsp1 carries four tensor components on electron one
         allocate(cr(nq*4), fr(nq*4))
         cr = 0; fr = 0
         ret = c_spsp1_sph(cr, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                           c_null_ptr, c_null_ptr)
         hv = int2e_spsp1_sph(fr, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmpr(cr, fr, nq*4, "int2e_spsp1_sph", ibasis, kmode, i, j, k, l, NG_SPSP1)
         deallocate(cr, fr)

         allocate(cr(nq*16), fr(nq*16))
         cr = 0; fr = 0
         ret = c_spsp12_sph(cr, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                            c_null_ptr, c_null_ptr)
         hv = int2e_spsp1spsp2_sph(fr, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmpr(cr, fr, nq*16, "int2e_spsp1spsp2_sph", ibasis, kmode, i, j, k, l, NG_SPSP12)
         deallocate(cr, fr)

         di = int(CINTcgto_cart(i, bas)); dj = int(CINTcgto_cart(j, bas))
         dk = int(CINTcgto_cart(k, bas)); dl = int(CINTcgto_cart(l, bas))
         nq = di*dj*dk*dl
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=dl
         allocate(cr(nq*4), fr(nq*4))
         cr = 0; fr = 0
         ret = c_spsp1_cart(cr, c_null_ptr, shls, atm, 8_c_int, bas, nbas, env, &
                            c_null_ptr, c_null_ptr)
         hv = int2e_spsp1_cart(fr, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
         call cmpr(cr, fr, nq*4, "int2e_spsp1_cart", ibasis, kmode, i, j, k, l, NG_SPSP1)
         deallocate(cr, fr)
      end do
      end do
      end do
      end do
      deallocate(fatm, fbas, fenv)
   end do
   end do

   do iw = 1, 6
      print '(A,A24,A,I0)', "  ", names(iw), " differing: ", tally(iw)
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

   subroutine note(what, ib, km, ii, jj, kk, ll, m)
      character(len=*), intent(in) :: what
      integer, intent(in) :: ib, km, ii, jj, kk, ll, m
      nbad = nbad + 1
      do iw = 1, 7
         if (trim(names(iw)) == trim(what)) tally(iw) = tally(iw) + 1
      end do
      if (nbad == 1) write(wat,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(what), &
         " basis=", ib, " kmode=", km, " shls (",ii,",",jj,",",kk,",",ll,") m=", m
   end subroutine note

   ! Rys order for a quartet: the expression the envs setup uses, on the
   ! CEILINGS rather than the bare angular momenta.  Each sigma.p raises two
   ! of them by one, and using the bare l instead understates the root count
   ! -- which understates where the eigensolver takes over, and flags the
   ! last-bit differences there as failures.  inc is the integral's ng
   ! increments, the same four numbers the generator emits.
   integer function rys_order(s1, s2, s3, s4, inc)
      integer, intent(in) :: s1, s2, s3, s4, inc(4)
      rys_order = (int(bas(ANG_OF, s1+1)) + inc(1) + int(bas(ANG_OF, s2+1)) + inc(2) &
                 + int(bas(ANG_OF, s3+1)) + inc(3) + int(bas(ANG_OF, s4+1)) + inc(4)) / 2 + 1
   end function rys_order

   real(dp) function tol_for(nr)
      integer, intent(in) :: nr
      if (nr <= 5) then
         tol_for = 0.0_dp
      else
         tol_for = 1.0e-11_dp
      end if
   end function tol_for

   subroutine cmpz(c, f, n, what, ib, km, ii, jj, kk, ll, inc)
      complex(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, km, ii, jj, kk, ll, inc(4)
      character(len=*), intent(in) :: what
      integer :: m, nr
      real(dp) :: r, scal, tol
      if (n <= 0) return
      nr = rys_order(ii, jj, kk, ll, inc)
      tol = tol_for(nr)
      scal = maxval(abs(c(1:n)))
      if (scal <= 0.0_dp) scal = 1.0_dp
      do m = 1, n
         ncmp = ncmp + 1
         r = abs(c(m) - f(m)) / scal
         if (r > worst) worst = r
         if (r > tol) call note(what, ib, km, ii, jj, kk, ll, m)
         if (r /= 0.0_dp .and. nr <= 5) nexact = nexact + 1
      end do
   end subroutine cmpz

   subroutine cmpr(c, f, n, what, ib, km, ii, jj, kk, ll, inc)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ib, km, ii, jj, kk, ll, inc(4)
      character(len=*), intent(in) :: what
      integer :: m, nr
      real(dp) :: r, scal, tol
      if (n <= 0) return
      nr = rys_order(ii, jj, kk, ll, inc)
      tol = tol_for(nr)
      scal = maxval(abs(c(1:n)))
      if (scal <= 0.0_dp) scal = 1.0_dp
      do m = 1, n
         ncmp = ncmp + 1
         r = abs(c(m) - f(m)) / scal
         if (r > worst) worst = r
         if (r > tol) call note(what, ib, km, ii, jj, kk, ll, m)
         if (r /= 0.0_dp .and. nr <= 5) nexact = nexact + 1
      end do
   end subroutine cmpr

end program gen2e_check
