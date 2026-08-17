!
! D10: the generated three- and two-centre catalogue against the C.
!
! Seventeen integrals in their real forms -- Cartesian and spherical -- from
! the same auto_intor.cl descriptions the C is generated from.  int3c2e and
! int2c2e are the controls: cint_3c2e.f90 already ships hand-written ones that
! int3c2e_check knows are bit-identical, so an emitter regression shows there
! first.
!
! THE COMPONENT COUNT OF EACH BLOCK IS READ OFF THE GENERATED ng ARRAY, not
! counted by eye.  int3c2e_spsp1 has four components rather than one, because
! the sigma.p pair makes ncomp_e1 four while dividing the tensor slot by the
! same factor -- and guessing it wrong overruns the buffer, which is a
! segfault if you are lucky and silence if you are not.  It was a segfault.
!
! int2c2e has no _spinor form here for the same reason it has none in the C:
! CINT2c2e_spinor_drv is an unimplemented stub upstream, and the C's own
! generator emits a function that prints and returns.
!
! Same bar as gen2e_check: bit-identical up to five Rys roots, 1e-11 above,
! with the root count from each integral's own ng increments.
!
program gen3c_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, ANG_OF, &
                                CINTcgto_cart, CINTcgto_spheric
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_gen_int3c2e
   use cint_3c2e, only: int3c2e_cart, int3c2e_sph, int2c2e_cart, int2c2e_sph
   implicit none

   interface
      function c_int3c2e_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ip1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ip1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ip1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ip2_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ip2_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ip2_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ip2_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_pvp1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_pvp1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_pvp1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_pvp1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_pvxp1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_pvxp1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_pvxp1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_pvxp1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_spsp1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_spsp1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_spsp1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_spsp1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ipspsp1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ipspsp1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ipspsp1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ipspsp1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_spsp1ip2_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_spsp1ip2_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_spsp1ip2_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_spsp1ip2_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ipip1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ipip1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ipip1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ipip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ipip2_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ipip2_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ipip2_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ipip2_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ipvip1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ipvip1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ipvip1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ipvip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ip1ip2_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ip1ip2_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c2e_ip1ip2_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_ip1ip2_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_ip1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_ip1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_ip1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_ip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_ip2_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_ip2_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_ip2_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_ip2_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_ipip1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_ipip1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_ipip1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_ipip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_ip1ip2_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_ip1ip2_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int2c2e_ip1ip2_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_ip1ip2_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
   end interface

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, shls(4), ret
   real(c_double), allocatable :: cb(:), fb(:)
   integer :: ib, di, dj, dk, nq, dims(0:3), fshls(0:3), sph, stride
   integer :: ncmp, nbad, nexact
   real(dp) :: worst
   character(len=80) :: wat
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv

   integer, parameter :: NG_INT3C2E(4) = [0,0,0,0]
   integer, parameter :: NG_INT3C2E_IP1(4) = [1,0,0,0]
   integer, parameter :: NG_INT3C2E_IP2(4) = [0,0,1,0]
   integer, parameter :: NG_INT3C2E_PVP1(4) = [1,1,0,0]
   integer, parameter :: NG_INT3C2E_PVXP1(4) = [1,1,0,0]
   integer, parameter :: NG_INT3C2E_SPSP1(4) = [1,1,0,0]
   integer, parameter :: NG_INT3C2E_IPSPSP1(4) = [2,1,0,0]
   integer, parameter :: NG_INT3C2E_SPSP1IP2(4) = [1,1,1,0]
   integer, parameter :: NG_INT3C2E_IPIP1(4) = [2,0,0,0]
   integer, parameter :: NG_INT3C2E_IPIP2(4) = [0,0,2,0]
   integer, parameter :: NG_INT3C2E_IPVIP1(4) = [1,1,0,0]
   integer, parameter :: NG_INT3C2E_IP1IP2(4) = [1,0,1,0]
   integer, parameter :: NG_INT2C2E(4) = [0,0,0,0]
   integer, parameter :: NG_INT2C2E_IP1(4) = [1,0,0,0]
   integer, parameter :: NG_INT2C2E_IP2(4) = [0,0,1,0]
   integer, parameter :: NG_INT2C2E_IPIP1(4) = [2,0,0,0]
   integer, parameter :: NG_INT2C2E_IP1IP2(4) = [1,0,1,0]

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ncmp = 0; nbad = 0; nexact = 0; worst = 0.0_dp; wat = "(none)"

   do ib = 1, n_reference_basis
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ib, bas, env, nbas, off)
      allocate(fatm(0:47), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1))
      do i = 1, 8
         fatm(6*(i-1):6*i-1) = int(atm(:, i))
      end do
      do i = 1, nbas
         fbas(8*(i-1):8*i-1) = int(bas(:, i))
      end do
      fenv = real(env, dp)

      stride = max(1, int(nbas)/5)
      do sph = 0, 1
      do i = 0, nbas-1, stride
      do j = 0, nbas-1, stride
      do k = 0, nbas-1, stride
         if (sph == 1) then
            di = int(CINTcgto_spheric(i,bas)); dj = int(CINTcgto_spheric(j,bas))
            dk = int(CINTcgto_spheric(k,bas))
         else
            di = int(CINTcgto_cart(i,bas)); dj = int(CINTcgto_cart(j,bas))
            dk = int(CINTcgto_cart(k,bas))
         end if
         shls(1)=i; shls(2)=j; shls(3)=k; shls(4)=0
         fshls(0)=i; fshls(1)=j; fshls(2)=k; fshls(3)=0

         nq = di*dj*dk*1
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e", ib, sph, i, j, k, NG_INT3C2E)
         deallocate(cb, fb)

         nq = di*dj*dk*3
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_ip1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_ip1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ip1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_ip1", ib, sph, i, j, k, NG_INT3C2E_IP1)
         deallocate(cb, fb)

         nq = di*dj*dk*3
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_ip2_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ip2_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_ip2_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ip2_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_ip2", ib, sph, i, j, k, NG_INT3C2E_IP2)
         deallocate(cb, fb)

         nq = di*dj*dk*1
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_pvp1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_pvp1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_pvp1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_pvp1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_pvp1", ib, sph, i, j, k, NG_INT3C2E_PVP1)
         deallocate(cb, fb)

         nq = di*dj*dk*3
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_pvxp1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_pvxp1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_pvxp1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_pvxp1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_pvxp1", ib, sph, i, j, k, NG_INT3C2E_PVXP1)
         deallocate(cb, fb)

         nq = di*dj*dk*4
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_spsp1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_spsp1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_spsp1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_spsp1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_spsp1", ib, sph, i, j, k, NG_INT3C2E_SPSP1)
         deallocate(cb, fb)

         nq = di*dj*dk*12
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_ipspsp1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ipspsp1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_ipspsp1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ipspsp1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_ipspsp1", ib, sph, i, j, k, NG_INT3C2E_IPSPSP1)
         deallocate(cb, fb)

         nq = di*dj*dk*12
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_spsp1ip2_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_spsp1ip2_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_spsp1ip2_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_spsp1ip2_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_spsp1ip2", ib, sph, i, j, k, NG_INT3C2E_SPSP1IP2)
         deallocate(cb, fb)

         nq = di*dj*dk*9
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_ipip1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ipip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_ipip1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ipip1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_ipip1", ib, sph, i, j, k, NG_INT3C2E_IPIP1)
         deallocate(cb, fb)

         nq = di*dj*dk*9
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_ipip2_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ipip2_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_ipip2_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ipip2_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_ipip2", ib, sph, i, j, k, NG_INT3C2E_IPIP2)
         deallocate(cb, fb)

         nq = di*dj*dk*9
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_ipvip1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ipvip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_ipvip1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ipvip1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_ipvip1", ib, sph, i, j, k, NG_INT3C2E_IPVIP1)
         deallocate(cb, fb)

         nq = di*dj*dk*9
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c2e_ip1ip2_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ip1ip2_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c2e_ip1ip2_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c2e_ip1ip2_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c2e_ip1ip2", ib, sph, i, j, k, NG_INT3C2E_IP1IP2)
         deallocate(cb, fb)

         nq = di*dj*1*1
         dims(0)=di; dims(1)=dj; dims(2)=1; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int2c2e_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int2c2e_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int2c2e", ib, sph, i, j, 0, NG_INT2C2E)
         deallocate(cb, fb)

         nq = di*dj*1*3
         dims(0)=di; dims(1)=dj; dims(2)=1; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int2c2e_ip1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_ip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int2c2e_ip1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_ip1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int2c2e_ip1", ib, sph, i, j, 0, NG_INT2C2E_IP1)
         deallocate(cb, fb)

         nq = di*dj*1*3
         dims(0)=di; dims(1)=dj; dims(2)=1; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int2c2e_ip2_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_ip2_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int2c2e_ip2_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_ip2_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int2c2e_ip2", ib, sph, i, j, 0, NG_INT2C2E_IP2)
         deallocate(cb, fb)

         nq = di*dj*1*9
         dims(0)=di; dims(1)=dj; dims(2)=1; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int2c2e_ipip1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_ipip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int2c2e_ipip1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_ipip1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int2c2e_ipip1", ib, sph, i, j, 0, NG_INT2C2E_IPIP1)
         deallocate(cb, fb)

         nq = di*dj*1*9
         dims(0)=di; dims(1)=dj; dims(2)=1; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int2c2e_ip1ip2_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_ip1ip2_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int2c2e_ip1ip2_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int2c2e_ip1ip2_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int2c2e_ip1ip2", ib, sph, i, j, 0, NG_INT2C2E_IP1IP2)
         deallocate(cb, fb)

      end do
      end do
      end do
      end do
      deallocate(fatm, fbas, fenv)
   end do

   print '(A,I0)',     "  values compared      : ", ncmp
   print '(A,ES10.2)', "  worst diff / block   : ", worst
   print '(A,I0)',     "  inexact at <=5 roots : ", nexact
   print '(A,I0)',     "  over tolerance       : ", nbad
   if (nbad > 0 .or. nexact > 0) then
      print '(A,A)',   "  worst at             : ", trim(wat)
      stop 1
   end if
   print '(A)', "  RESULT: PASS"

contains

   integer function rys_order(s1, s2, s3, ic)
      integer, intent(in) :: s1, s2, s3, ic(4)
      rys_order = (int(bas(ANG_OF,s1+1)) + ic(1) + int(bas(ANG_OF,s2+1)) + ic(2) &
                 + int(bas(ANG_OF,s3+1)) + ic(3) + ic(4)) / 2 + 1
   end function rys_order

   subroutine cmp(c, f, n, what, ibb, sp, ii, jj, kk, ic)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ibb, sp, ii, jj, kk, ic(4)
      character(len=*), intent(in) :: what
      integer :: m, nr
      real(dp) :: r, scal, tol
      if (n <= 0) return
      nr = rys_order(ii, jj, kk, ic)
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
            if (nbad == 1) write(wat,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(what), &
               " basis=", ibb, " sph=", sp, " shls (",ii,",",jj,",",kk,") m=", m, " roots=", nr
         end if
         if (r /= 0.0_dp .and. nr <= 5) nexact = nexact + 1
      end do
   end subroutine cmp

end program gen3c_check
