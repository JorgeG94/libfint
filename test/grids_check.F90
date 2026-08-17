!
! D10: int1e_grids against the C.
!
! The electrostatic-potential integrals: <i| 1/|r-r_g| |j> at a list of points
! rather than at nuclei.  This is the one family whose g array has a shape of
! its own -- the grid index is innermost and blocked -- so it gets its own
! check rather than a row in the catalogue sweep.
!
! The grid is deliberately awkward.  GRID_BLKSIZE is 104, and everything
! interesting in this code happens at a block boundary: the last block is
! short, `bgrids` differs from GRID_BLKSIZE there, and the contraction and the
! output copy both index by it.  So the counts swept are 1, 103, 104, 105 and
! 209 -- one, one short of a block, exactly a block, one past, and two blocks
! plus one.  A sweep that used a round multiple would exercise none of that.
!
! The points themselves are put where they can catch sign and ordering errors:
! on top of a nucleus is skipped (the integral diverges), but near one, far
! away, and off-axis are all included.
!
! The four generated derivatives -- ip, ipvip, spvsp, ipip -- ride along, so
! the grids gout emitter is checked at the same time as the runtime it sits on.
!
! The bar is bit-identical.
!
program grids_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, &
                                CINTcgto_cart, CINTcgto_spheric
   use cint_test_systems
   use cint_const, only: dp
   use cint_workspace, only: cint_ws
   use cint_1e_grids, only: int1e_grids_cart, int1e_grids_sph, GRID_BLKSIZE
   use cint_gen_int1e_grids1
   implicit none

   interface
      function c_grids_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_grids_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_ip_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_ip_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_ip_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_ip_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_ipvip_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_ipvip_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_ipvip_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_ipvip_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_spvsp_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_spvsp_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_spvsp_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_spvsp_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_ipip_cart(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_ipip_cart') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
      function c_int1e_grids_ipip_sph(o,d,s,a,na,b,nb,e,p,c) bind(C,name='int1e_grids_ipip_sph') result(r)
         import :: c_int, c_double, c_ptr
         real(c_double)::o(*); type(c_ptr),value::d,p,c
         integer(c_int)::s(*),a(*),b(*); integer(c_int),value::na,nb
         real(c_double)::e(*); integer(c_int)::r
      end function
   end interface

   ! libcint_interface has no name for it; the C's is include/cint.h.in:44,
   ! and this is the 1-based Fortran slot for that 0-based C offset.
   integer, parameter :: PTR_GRIDS_F = 12 + 1

   integer, parameter :: NG_CASES = 5
   integer, parameter :: GRID_COUNTS(NG_CASES) = [1, 103, 104, 105, 209]

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, shls(4), ret
   real(c_double), allocatable :: cb(:), fb(:)
   integer :: ib, ic, ng, di, dj, nq, dims(0:3), fshls(0:3), sph, stride, g
   integer :: ncmp, nbad, gridptr
   character(len=80) :: wat
   type(cint_ws) :: ws
   integer, allocatable :: fatm(:), fbas(:)
   real(dp), allocatable :: fenv(:)
   logical :: hv
   real(dp) :: t

   allocate(atm(ATM_SLOTS, 8), bas(BAS_SLOTS, 200), env(60000))
   ncmp = 0; nbad = 0; wat = "(none)"

   do ib = 1, 3
   do ic = 1, NG_CASES
      ng = GRID_COUNTS(ic)
      atm = 0; bas = 0; env = 0.0_c_double
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ib, bas, env, nbas, off)

      ! The grid lives in env, and PTR_GRIDS holds its 0-based offset -- so
      ! the value stored is `off`, while PTR_GRIDS_F is a 1-based slot index
      ! into the Fortran array.  Getting those two the wrong way round is the
      ! D0 bug in a new place.
      !
      ! `off` is only the free offset because cint_test_systems now advances
      ! it past the basis block.  It did not, and the grid landed on top of
      ! the first shell's exponents -- which made both libcint and the port
      ! return NaN, in perfect agreement, for a reason that had nothing to do
      ! with either.
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
      do sph = 0, 1
      do i = 0, nbas-1, stride
      do j = 0, nbas-1, stride
         if (sph == 1) then
            di = int(CINTcgto_spheric(i,bas)); dj = int(CINTcgto_spheric(j,bas))
         else
            di = int(CINTcgto_cart(i,bas)); dj = int(CINTcgto_cart(j,bas))
         end if
         nq = di*dj*ng
         ! shls(3) and shls(4) are the grid range, not shells
         shls(1)=i; shls(2)=j; shls(3)=0; shls(4)=ng
         fshls(0)=i; fshls(1)=j; fshls(2)=0; fshls(3)=ng
         dims(0)=di; dims(1)=dj; dims(2)=ng; dims(3)=1
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_grids_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_grids_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, ib, ng, sph, i, j)
         deallocate(cb, fb)

         nq = di*dj*ng*3
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int1e_grids_ip_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_ip_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int1e_grids_ip_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_ip_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, ib, ng, sph, i, j)
         deallocate(cb, fb)

         nq = di*dj*ng*9
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int1e_grids_ipvip_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_ipvip_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int1e_grids_ipvip_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_ipvip_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, ib, ng, sph, i, j)
         deallocate(cb, fb)

         nq = di*dj*ng*4
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int1e_grids_spvsp_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_spvsp_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int1e_grids_spvsp_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_spvsp_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, ib, ng, sph, i, j)
         deallocate(cb, fb)

         nq = di*dj*ng*9
         allocate(cb(nq), fb(nq)); cb=0; fb=0
         if (sph == 1) then
            ret = c_int1e_grids_ipip_sph(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_ipip_sph(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         else
            ret = c_int1e_grids_ipip_cart(cb,c_null_ptr,shls,atm,8_c_int,bas,nbas,env,c_null_ptr,c_null_ptr)
            hv = int1e_grids_ipip_cart(fb,dims,fshls,fatm,8,fbas,int(nbas),fenv,ws)
         end if
         call cmp(cb, fb, nq, ib, ng, sph, i, j)
         deallocate(cb, fb)
      end do
      end do
      end do
      deallocate(fatm, fbas, fenv)
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

   subroutine cmp(c, f, n, ibb, ngg, sp, ii, jj)
      real(dp), intent(in) :: c(:), f(:)
      integer, intent(in) :: n, ibb, ngg, sp, ii, jj
      integer :: m
      do m = 1, n
         ncmp = ncmp + 1
         if (c(m) /= f(m)) then
            nbad = nbad + 1
            if (nbad == 1) then
               write(wat,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') &
                  "basis=", ibb, " ngrids=", ngg, " sph=", sp, &
                  " shls (", ii, ",", jj, ") m=", m
               print '(A,I0,A,2ES24.15)', "   n=", n, "  C, port = ", c(m), f(m)
            end if
         end if
      end do
   end subroutine cmp

end program grids_check
