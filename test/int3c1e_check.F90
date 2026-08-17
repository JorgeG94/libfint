!
! D10: the three-centre one-electron integrals against the C.
!
! int3c1e (a plain three-centre overlap) and int3c1e_rinv (the same with a
! 1/r kernel between the pair and the auxiliary shell, so a Rys sum).
!
! Worth its own check rather than a row in the catalogue sweep because its g
! array is built along j rather than i, and the j stride changes partway
! through the recursion -- from li+1 in the compressed (i,j) plane to
! g_stride_j in the full (i,j,k) block.  Nothing else in the port does that.
!
! The three generated descriptions -- p2, iprinv, ip1 -- ride along, so the
! 3c1e gout emitter is checked with the runtime it sits on.
!
! The bar is bit-identical.
!
program int3c1e_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, &
                                CINTcgto_cart, CINTcgto_spheric
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_3c1e, only: int3c1e_cart, int3c1e_sph, &
                        int3c1e_rinv_cart, int3c1e_rinv_sph
   use cint_gen_int3c1e
   implicit none
   interface
      function c_o_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_o_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_r_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_rinv_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_r_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_rinv_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c1e_p2_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_p2_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c1e_p2_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_p2_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c1e_iprinv_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_iprinv_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c1e_iprinv_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_iprinv_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c1e_ip1_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_ip1_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int3c1e_ip1_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c1e_ip1_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
   end interface

   integer, parameter :: PTR_RINV_ORIG_F = 4 + 1

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, shls(4), ret
   real(c_double), allocatable :: cb(:), fb(:)
   integer :: ib, di, dj, dk, nq, dims(0:3), fshls(0:3), sph, stride
   integer :: ncmp, nbad
   character(len=80) :: wat
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(20000))
   ncmp = 0; nbad = 0; wat = "(none)"

   do ib = 1, n_reference_basis
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ib, bas, env, nbas, off)
      ! somewhere off every nucleus, so rinv is not evaluated at a singularity
      env(PTR_RINV_ORIG_F+0) = 0.21_c_double
      env(PTR_RINV_ORIG_F+1) = -0.34_c_double
      env(PTR_RINV_ORIG_F+2) = 0.55_c_double

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
         nq = di*dj*dk
         shls(1)=i; shls(2)=j; shls(3)=k; shls(4)=0
         fshls(0)=i; fshls(1)=j; fshls(2)=k; fshls(3)=0
         dims(0)=di; dims(1)=dj; dims(2)=dk; dims(3)=1

         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_o_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_o_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c1e", ib, sph, i, j, k)
         deallocate(cb, fb)

         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_r_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_rinv_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_r_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_rinv_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, "int3c1e_rinv", ib, sph, i, j, k)
         deallocate(cb, fb)

         allocate(cb(nq*1), fb(nq*1)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c1e_p2_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_p2_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c1e_p2_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_p2_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq*1, "int3c1e_p2", ib, sph, i, j, k)
         deallocate(cb, fb)

         allocate(cb(nq*3), fb(nq*3)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c1e_iprinv_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_iprinv_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c1e_iprinv_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_iprinv_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq*3, "int3c1e_iprinv", ib, sph, i, j, k)
         deallocate(cb, fb)

         allocate(cb(nq*3), fb(nq*3)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int3c1e_ip1_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_ip1_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int3c1e_ip1_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int3c1e_ip1_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq*3, "int3c1e_ip1", ib, sph, i, j, k)
         deallocate(cb, fb)
      end do
      end do
      end do
      end do
      deallocate(fatm, fbas, fenv)
   end do

   print '(A,I0)', "  values compared : ", ncmp
   print '(A,I0)', "  differing       : ", nbad
   if (nbad > 0) then
      print '(A,A)', "  first at        : ", trim(wat)
      stop 1
   end if
   print '(A)', "  RESULT: PASS (bit-identical)"

contains

   subroutine cmp(c, f, n, what, ibb, sp, ii, jj, kk)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ibb, sp, ii, jj, kk
      character(len=*), intent(in) :: what
      integer :: m
      do m = 1, n
         ncmp = ncmp + 1
         if (c(m) /= f(m)) then
            nbad = nbad + 1
            if (nbad == 1) write(wat,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(what), &
               " basis=", ibb, " sph=", sp, " shls (",ii,",",jj,",",kk,") m=", m
         end if
      end do
   end subroutine cmp

end program int3c1e_check
