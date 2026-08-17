!
! D9: the three-centre spinor integrals against the C.
!
! int3c2e_spinor and its ssc twin, over the shell triples of all five
! reference bases in all three kappa modes.
!
! ssc -- "spheric-spinor-cartesian" -- leaves the auxiliary index Cartesian,
! which is the form a density-fitted relativistic code wants because the
! fitting basis has no spin to carry.  It is worth testing separately and not
! only for coverage: the ssc path skips the inner spherical transform
! entirely, so it is the one case where getting that transform's "did it
! fire?" sentinel backwards produces right answers.  Which is exactly what
! happened -- RESULT_IN_GCART is zero, a bare `loc = 0` said the opposite of
! what it looked like, and only the non-ssc half noticed.
!
! The bar is bit-identical: three-centre integrals over these bases stay
! inside the closed-form Rys roots.
!
program spinor3c_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, ANG_OF, KAPPA_OF, &
                                CINTcgto_spinor, CINTcgto_spheric, CINTcgto_cart
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_3c2e_spinor, only: int3c2e_spinor, int3c2e_spinor_ssc, int2c2e_spinor
   implicit none
   interface
      function c3(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c3s(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int3c2e_spinor_ssc') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c2c(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int2c2e_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
   end interface
   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, shls(4), ret
   complex(c_double_complex), allocatable :: cb(:), fb(:)
   integer :: ib,km,di,dj,dk,nq,dims(0:3),fshls(0:3),m,lang,stride,ncmp,nbad,nc2,nb2
   integer :: nc3, nb3
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   type(cint_ws) :: ws
   logical :: hv
   allocate(atm(ATM_SLOTS,8), bas(BAS_SLOTS,200), env(20000))
   ncmp=0; nbad=0; nc2=0; nb2=0; nc3=0; nb3=0
   do ib=1,n_reference_basis
   do km=0,2
      atm=0; bas=0; env=0
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ib, bas, env, nbas, off)
      do i=1,nbas
         lang=int(bas(ANG_OF,i))
         select case(km)
         case(1); bas(KAPPA_OF,i)=int(-lang-1,c_int)
         case(2); if(lang>0) bas(KAPPA_OF,i)=int(lang,c_int)
         end select
      end do
      allocate(fatm(0:47), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1))
      do i=1,8; fatm(6*(i-1):6*i-1)=int(atm(:,i)); end do
      do i=1,nbas; fbas(8*(i-1):8*i-1)=int(bas(:,i)); end do
      fenv = real(env, dp)
      stride = max(1,int(nbas)/6)
      do i=0,nbas-1,stride
      do j=0,nbas-1,stride
      do k=0,nbas-1,stride
         di=int(CINTcgto_spinor(i,bas)); dj=int(CINTcgto_spinor(j,bas))
         dk=int(CINTcgto_spheric(k,bas))
         if (di<=0 .or. dj<=0 .or. dk<=0) cycle
         nq=di*dj*dk
         shls=[i,j,k,0]; fshls=shls; dims=[di,dj,dk,1]
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         ret = c3(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int3c2e_spinor(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         do m=1,nq; ncmp=ncmp+1; if (cb(m)/=fb(m)) nbad=nbad+1; end do
         deallocate(cb,fb)
         ! ssc: the auxiliary index stays Cartesian
         dk = int(CINTcgto_cart(k, bas))
         nq=di*dj*dk
         dims=[di,dj,dk,1]
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         ret = c3s(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int3c2e_spinor_ssc(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         do m=1,nq; nc2=nc2+1; if (cb(m)/=fb(m)) nb2=nb2+1; end do
         deallocate(cb,fb)
      end do; end do; end do
      ! int2c2e_spinor.  Two centres, so its own pair loop.  CINT2c2e_spinor_drv
      ! is the only spinor two-centre driver upstream implements: the four
      ! derivative forms print "&c2s_sf_1e_spinor not implemented" and return
      ! zero, so there is nothing there to compare against.
      do i=0,nbas-1,stride
      do k=0,nbas-1,stride
         di=int(CINTcgto_spinor(i,bas)); dk=int(CINTcgto_spinor(k,bas))
         if (di<=0 .or. dk<=0) cycle
         nq=di*dk
         shls=[i,k,0,0]; fshls=shls; dims=[di,dk,1,1]
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         ret = c2c(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int2c2e_spinor(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         do m=1,nq; nc3=nc3+1; if (cb(m)/=fb(m)) nb3=nb3+1; end do
         deallocate(cb,fb)
      end do; end do
      deallocate(fatm,fbas,fenv)
   end do; end do
   print '(A,I0,A,I0)', "  int3c2e_spinor      values=", ncmp, "  differing=", nbad
   print '(A,I0,A,I0)', "  int3c2e_spinor_ssc  values=", nc2,  "  differing=", nb2
   print '(A,I0,A,I0)', "  int2c2e_spinor      values=", nc3,  "  differing=", nb3
   if (nbad>0 .or. nb2>0 .or. nb3>0) stop 1
   print '(A)', "  RESULT: PASS (bit-identical)"
end program spinor3c_check
