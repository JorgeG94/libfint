!
! The one-electron grids spinor integrals against the C.
!
! Five entry points -- int1e_grids_spinor and the four generated derivatives --
! over the awkward grid counts grids_check already established are the ones
! that matter: 1, 103, 104, 105 and 209, which straddle GRID_BLKSIZE in every
! direction that changes `bgrids`.
!
! Four of the five use the spin-free transform and one, spvsp, the
! spin-included one, so both halves of apply_c2s_spinor_1e_grids are exercised.
! All three kappa modes are swept, because the spinor length depends on kappa
! and not only on l, and the transform indexes its coefficient table by both.
!
! The bar is bit-identical.
!
program spinor_grids_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, ANG_OF, &
                                KAPPA_OF, CINTcgto_spinor
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_1e_grids, only: GRID_BLKSIZE
   use cint_1e_grids_spinor, only: int1e_grids_spinor
   use cint_gen_int1e_grids1
   implicit none
   interface
      function c_int1e_grids(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_ip(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_ip_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_ipvip(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_ipvip_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_spvsp(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_spvsp_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_ipip(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_ipip_spinor') result(r)
         import :: c_int, c_double, c_double_complex, c_ptr
         complex(c_double_complex)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
   end interface

   ! include/cint.h.in:44 as a 1-based Fortran slot
   integer, parameter :: PTR_GRIDS_F = 12 + 1
   integer, parameter :: NG_CASES = 5
   integer, parameter :: GRID_COUNTS(NG_CASES) = [1, 103, 104, 105, 209]

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, shls(4), ret
   complex(c_double_complex), allocatable :: cb(:), fb(:)
   integer :: ib, ic, ng, di, dj, nq, dims(0:3), fshls(0:3), stride, g, km, lang
   integer :: ncmp, nbad, gridptr
   character(len=90) :: wat
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv
   real(dp) :: t

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(60000))
   ncmp = 0; nbad = 0; wat = "(none)"

   do ib = 1, 3
   do ic = 1, NG_CASES
   do km = 0, 2
      ng = GRID_COUNTS(ic)
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ib, bas, env, nbas, off)
      do i = 1, nbas
         lang = int(bas(ANG_OF,i))
         select case(km)
         case(1); bas(KAPPA_OF,i) = int(-lang-1, c_int)
         case(2); if (lang > 0) bas(KAPPA_OF,i) = int(lang, c_int)
         end select
      end do

      gridptr = int(off)
      env(PTR_GRIDS_F) = real(gridptr, c_double)
      do g = 0, ng-1
         t = real(g, dp)
         env(off+1) = 0.31_dp * cos(t*0.7_dp) + 0.05_dp * t
         env(off+2) = 0.47_dp * sin(t*1.3_dp) - 0.03_dp * t
         env(off+3) = 0.93_dp + 0.11_dp * cos(t*2.1_dp)
         off = off + 3
      end do

      allocate(fatm(0:47), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1))
      do i = 1, 8
         fatm(6*(i-1):6*i-1) = int(atm(:, i))
      end do
      do i = 1, nbas
         fbas(8*(i-1):8*i-1) = int(bas(:, i))
      end do
      fenv = real(env, dp)

      stride = max(1, int(nbas)/4)
      do i = 0, nbas-1, stride
      do j = 0, nbas-1, stride
         di = int(CINTcgto_spinor(i, bas)); dj = int(CINTcgto_spinor(j, bas))
         if (di <= 0 .or. dj <= 0) cycle
         ! shls(3) and shls(4) are the grid range, not shells
         shls(1)=i; shls(2)=j; shls(3)=0; shls(4)=ng
         fshls(0)=i; fshls(1)=j; fshls(2)=0; fshls(3)=ng
         dims(0)=di; dims(1)=dj; dims(2)=ng; dims(3)=1

         nq = di*dj*ng*1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         ret = c_int1e_grids(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int1e_grids_spinor(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         call cmpd(cb, fb, nq, ib, ng, km, i, j, "int1e_grids")

         nq = di*dj*ng*3
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         ret = c_int1e_grids_ip(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int1e_grids_ip_spinor(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         call cmpd(cb, fb, nq, ib, ng, km, i, j, "int1e_grids_ip")

         nq = di*dj*ng*9
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         ret = c_int1e_grids_ipvip(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int1e_grids_ipvip_spinor(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         call cmpd(cb, fb, nq, ib, ng, km, i, j, "int1e_grids_ipvip")

         nq = di*dj*ng*4
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         ret = c_int1e_grids_spvsp(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int1e_grids_spvsp_spinor(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         call cmpd(cb, fb, nq, ib, ng, km, i, j, "int1e_grids_spvsp")

         nq = di*dj*ng*9
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         ret = c_int1e_grids_ipip(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
         hv = int1e_grids_ipip_spinor(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         call cmpd(cb, fb, nq, ib, ng, km, i, j, "int1e_grids_ipip")
      end do
      end do
      deallocate(fatm, fbas, fenv)
   end do
   end do
   end do

   print '(A,I0)', "  block size      : ", GRID_BLKSIZE
   print '(A,I0)', "  values compared : ", ncmp
   print '(A,I0)', "  differing       : ", nbad
   if (nbad > 0) then
      print '(A,A)', "  first at        : ", trim(wat)
      stop 1
   end if
   print '(A)', "  RESULT: PASS (bit-identical)"

contains

   ! compare, then release: the call sites allocate a fresh pair per integral
   ! and there are five of them per shell pair
   subroutine cmpd(c, f, n, ibb, ngg, kmm, ii, jj, what)
      complex(c_double_complex), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ibb, ngg, kmm, ii, jj
      character(len=*), intent(in) :: what
      call cmp(c, f, n, ibb, ngg, kmm, ii, jj, what)
      deallocate(cb, fb)
   end subroutine cmpd

   subroutine cmp(c, f, n, ibb, ngg, kmm, ii, jj, what)
      complex(c_double_complex), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ibb, ngg, kmm, ii, jj
      character(len=*), intent(in) :: what
      integer :: m
      do m = 1, n
         ncmp = ncmp + 1
         if (c(m) /= f(m)) then
            nbad = nbad + 1
            if (nbad == 1) then
               write(wat,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') trim(what), &
                  " basis=", ibb, " ngrids=", ngg, " kappa=", kmm, &
                  " shls (", ii, ",", jj, ") m=", m
               print '(A,I0,A,4ES22.13)', "   n=", n, "  C, port = ", &
                  real(c(m),dp), aimag(c(m)), real(f(m),dp), aimag(f(m))
            end if
         end if
      end do
   end subroutine cmp

end program spinor_grids_check
