!
! The F12 families against the C.
!
! int2e_yp and int2e_stg -- the Yukawa potential and the Slater-type geminal
! -- with their four derivative forms each.  Ten entry points, spherical
! only, which is all the C emits for this family.
!
! Two regimes, because the kernels branch on the geminal exponent.  With
! zeta > 0 the quadrature comes from CINTstg_roots, a two-dimensional
! Chebyshev interpolation over a 3.5-million-coefficient table; with zeta at
! or below zero both kernels fall back to the ordinary Rys roots and skip
! the weight rescaling, which is how the C spells "no geminal".  Only the
! first regime exercises the STG table at all.
!
! The bar is bit-identical.  The STG quadrature is table lookup and Clenshaw
! recurrence -- no eigensolver, so none of D3's DSTEMR deviation applies.
!
program f12_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, &
                                CINTcgto_spheric
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_gen_f12
   implicit none
   interface
      function c_yp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_yp_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_stg(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_stg_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_yp_ip1(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_yp_ip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_stg_ip1(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_stg_ip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_yp_ipip1(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_yp_ipip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_stg_ipip1(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_stg_ipip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_yp_ipvip1(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_yp_ipvip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_stg_ipvip1(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_stg_ipvip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_yp_ip1ip2(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_yp_ip1ip2_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_stg_ip1ip2(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2e_stg_ip1ip2_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
   end interface

   ! include/cint.h.in:40, as a 1-based Fortran slot
   integer, parameter :: PTR_F12_ZETA_F = 9 + 1

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, l, shls(4), ret
   real(dp), allocatable :: cb(:), fb(:)
   integer :: ib, iz, di, dj, dk, dl, nq, dims(0:3), fshls(0:3), m, stride
   integer :: ncmp, nbad
   character(len=90) :: wat
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   type(cint_ws) :: ws
   logical :: hv
   real(dp), parameter :: ZETAS(2) = [1.2_dp, 0.0_dp]

   allocate(atm(ATM_SLOTS,8), bas(BAS_SLOTS,200), env(20000))
   allocate(cb(400000), fb(400000))
   ncmp = 0; nbad = 0; wat = "(none)"

   do ib = 1, 3
   do iz = 1, 2
      atm=0; bas=0; env=0
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ib, bas, env, nbas, off)
      env(PTR_F12_ZETA_F) = ZETAS(iz)

      allocate(fatm(0:47), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1))
      do i=1,8; fatm(6*(i-1):6*i-1)=int(atm(:,i)); end do
      do i=1,nbas; fbas(8*(i-1):8*i-1)=int(bas(:,i)); end do
      fenv = real(env, dp)

      stride = max(1, int(nbas)/4)
      do i=0,nbas-1,stride
      do j=0,nbas-1,stride
      do k=0,nbas-1,stride
      do l=0,nbas-1,stride
         di=int(CINTcgto_spheric(i,bas)); dj=int(CINTcgto_spheric(j,bas))
         dk=int(CINTcgto_spheric(k,bas)); dl=int(CINTcgto_spheric(l,bas))
         nq = di*dj*dk*dl
         shls=[i,j,k,l]; fshls=shls; dims=[di,dj,dk,dl]

         call one(1, nq)
         call one(2, nq)
         call one(3, nq*3)
         call one(4, nq*3)
         call one(5, nq*9)
         call one(6, nq*9)
         call one(7, nq*9)
         call one(8, nq*9)
         call one(9, nq*9)
         call one(10, nq*9)
      end do; end do; end do; end do
      deallocate(fatm,fbas,fenv)
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

   subroutine one(which, n)
      integer, intent(in) :: which, n
      integer :: mm
      cb(1:n) = 0; fb(1:n) = 0
      select case (which)
      case (1)
         ret = c_yp(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_yp_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      case (2)
         ret = c_stg(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_stg_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      case (3)
         ret = c_yp_ip1(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_yp_ip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      case (4)
         ret = c_stg_ip1(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_stg_ip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      case (5)
         ret = c_yp_ipip1(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_yp_ipip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      case (6)
         ret = c_stg_ipip1(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_stg_ipip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      case (7)
         ret = c_yp_ipvip1(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_yp_ipvip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      case (8)
         ret = c_stg_ipvip1(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_stg_ipvip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      case (9)
         ret = c_yp_ip1ip2(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_yp_ip1ip2_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      case (10)
         ret = c_stg_ip1ip2(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2e_stg_ip1ip2_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
      end select
      do mm = 1, n
         ncmp = ncmp + 1
         if (cb(mm) /= fb(mm)) then
            nbad = nbad + 1
            if (nbad == 1) then
               write(wat,'(A,I0,A,F4.1,A,I0,A,4I4,A,I0)') "kernel=", which, &
                  " zeta=", ZETAS(iz), " basis=", ib, " shls", i, j, k, l, " m=", mm
               print '(A,2ES24.15)', "   C, port = ", cb(mm), fb(mm)
            end if
         end if
      end do
   end subroutine one

end program f12_check
