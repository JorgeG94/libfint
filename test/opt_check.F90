!
! The optimizer must not move a number.
!
! CINTOpt caches what the drivers would otherwise recompute for every shell
! quartet.  Three of its four caches -- the logged coefficient bounds, the
! non-zero contraction map, the Cartesian index map -- hold exactly the
! values the driver would have computed, so every integral has to come back
! bit-identical.  This checks that over all six arities the port has an
! optimizer for: 1e, 2e, 3c2e, 2c2e, 3c1e and 1e_grids.
!
! The fourth cache, the primitive-pair data, is deliberately not used by the
! drivers.  The C builds it with a shared angular bound -- li+ijkl_inc rather
! than the exact li_ceil, lj_ceil -- which is looser, keeps more primitives,
! and so gives the C's own optimised and unoptimised paths different last
! bits for any integral carrying derivatives on both sides.  Inheriting the
! cache without inheriting that would be worse than not caching it.
!
! The optimizer is also checked for the thing it must refuse: one built for
! a different arity, or a different ng, must not be used.  A silent mismatch
! there reads the wrong index map and returns plausible wrong numbers, which
! is exactly the failure the C has no guard against.
!
program opt_check
   use iso_c_binding
   use libcint_interface, only: ATM_SLOTS, BAS_SLOTS, PTR_ENV_START, &
                                CINTcgto_spheric
   use cint_test_systems
   use cint_const, only: dp
   use cint_envs,  only: cint_opt_t
   use cint_workspace, only: cint_ws
   use cint_opt,      only: cint_del_optimizer
   use cint_g1e,      only: cint_all_1e_optimizer
   use cint_g2e,      only: cint_all_2e_optimizer, cint_all_3c2e_optimizer, &
                            cint_all_2c2e_optimizer
   use cint_1e,       only: int1e_ovlp_sph
   use cint_2e,       only: int2e_sph
   use cint_3c2e,     only: int3c2e_sph, int2c2e_sph
   use cint_3c1e,     only: int3c1e_sph, cint_all_3c1e_optimizer
   use cint_1e_grids, only: int1e_grids_sph, cint_all_1e_grids_optimizer
   implicit none

   integer, parameter :: PTR_GRIDS_F = 12 + 1
   integer, parameter :: NG_OVLP(0:7)  = [0, 0, 0, 0, 0, 1, 1, 1]
   integer, parameter :: NG_2E(0:7)    = [0, 0, 0, 0, 0, 1, 1, 1]
   integer, parameter :: NG_3C2E(0:7)  = [0, 0, 0, 0, 0, 1, 1, 1]
   integer, parameter :: NG_2C2E(0:7)  = [0, 0, 0, 0, 0, 1, 1, 1]
   integer, parameter :: NG_3C1E(0:7)  = [0, 0, 0, 0, 0, 1, 0, 1]
   integer, parameter :: NG_GRIDS(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
   integer, parameter :: NGRID = 105

   integer(c_int), allocatable :: atm(:,:), bas(:,:)
   real(c_double), allocatable :: env(:)
   integer(c_int) :: nbas, off, i, j, k, l
   real(dp), allocatable :: b0(:), b1(:)
   integer :: ib, dims(0:3), nq, di, dj, dk, dl, m, stride, g, gridptr
   integer, allocatable :: fatm(:), fbas(:), fshls(:)
   real(dp), allocatable :: fenv(:)
   type(cint_ws) :: ws
   type(cint_opt_t), target :: o1e, o2e, o3c2e, o2c2e, o3c1e, ogr
   integer :: n(6), nb(6)
   real(dp) :: t
   logical :: hv
   character(len=12), parameter :: NAMES(6) = &
      [character(len=12) :: "int1e_ovlp","int2e","int3c2e","int2c2e","int3c1e","int1e_grids"]

   allocate(atm(ATM_SLOTS,8), bas(BAS_SLOTS,200), env(60000))
   allocate(b0(400000), b1(400000))
   n = 0; nb = 0

   do ib = 1, n_reference_basis
      atm=0; bas=0; env=0
      off = PTR_ENV_START
      call setup_c2h6_geometry(atm, env, off)
      call setup_reference_basis(ib, bas, env, nbas, off)

      gridptr = int(off)
      env(PTR_GRIDS_F) = real(gridptr, c_double)
      do g = 0, NGRID-1
         t = real(g, dp)
         env(off+1) = 0.31_dp * cos(t*0.7_dp) + 0.05_dp * t
         env(off+2) = 0.47_dp * sin(t*1.3_dp) - 0.03_dp * t
         env(off+3) = 0.93_dp + 0.11_dp * cos(t*2.1_dp)
         off = off + 3
      end do

      allocate(fatm(0:47), fbas(0:8*int(nbas)-1), fenv(0:size(env)-1), fshls(0:3))
      do i=1,8; fatm(6*(i-1):6*i-1)=int(atm(:,i)); end do
      do i=1,nbas; fbas(8*(i-1):8*i-1)=int(bas(:,i)); end do
      fenv = real(env, dp)

      call cint_all_1e_optimizer      (o1e,   NG_OVLP,  fatm, 8, fbas, int(nbas), fenv)
      call cint_all_2e_optimizer      (o2e,   NG_2E,    fatm, 8, fbas, int(nbas), fenv)
      call cint_all_3c2e_optimizer    (o3c2e, NG_3C2E,  fatm, 8, fbas, int(nbas), fenv)
      call cint_all_2c2e_optimizer    (o2c2e, NG_2C2E,  fatm, 8, fbas, int(nbas), fenv)
      call cint_all_3c1e_optimizer    (o3c1e, NG_3C1E,  fatm, 8, fbas, int(nbas), fenv)
      call cint_all_1e_grids_optimizer(ogr,   NG_GRIDS, fatm, 8, fbas, int(nbas), fenv)

      stride = max(1, int(nbas)/5)
      do i=0,nbas-1,stride
      do j=0,nbas-1,stride
         di=int(CINTcgto_spheric(i,bas)); dj=int(CINTcgto_spheric(j,bas))

         ! 1e
         fshls(0)=i; fshls(1)=j; fshls(2)=0; fshls(3)=0
         dims=[di,dj,1,1]
         call both(1, o1e, di*dj)

         ! 1e_grids -- shls(2), shls(3) are the grid range, not shells
         fshls(2)=0; fshls(3)=NGRID
         dims=[di,dj,NGRID,1]
         call both(6, ogr, di*dj*NGRID)

         ! 2c2e
         fshls(0)=i; fshls(1)=j; fshls(2)=0; fshls(3)=0
         dims=[di,dj,1,1]
         call both(4, o2c2e, di*dj)

         do k=0,nbas-1,stride
            dk=int(CINTcgto_spheric(k,bas))
            fshls(0)=i; fshls(1)=j; fshls(2)=k; fshls(3)=0
            dims=[di,dj,dk,1]
            call both(3, o3c2e, di*dj*dk)
            call both(5, o3c1e, di*dj*dk)
            do l=0,nbas-1,stride
               dl=int(CINTcgto_spheric(l,bas))
               fshls(3)=l
               dims=[di,dj,dk,dl]
               call both(2, o2e, di*dj*dk*dl)
            end do
         end do
      end do
      end do

      ! An optimizer for the wrong arity must be refused, not used.  2e shells
      ! through the 3c2e optimizer: same basis, same ng, wrong order.
      fshls=[0,0,0,0]
      di=int(CINTcgto_spheric(0,bas))
      dims=[di,di,di,di]
      nq = di**4
      b0=0; b1=0
      ws%opt => null()
      hv = int2e_sph(b0, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
      ws%opt => o3c2e
      hv = int2e_sph(b1, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
      ws%opt => null()
      do m=1,nq
         if (b0(m) /= b1(m)) nb(2) = nb(2) + 1
      end do

      call cint_del_optimizer(o1e);   call cint_del_optimizer(o2e)
      call cint_del_optimizer(o3c2e); call cint_del_optimizer(o2c2e)
      call cint_del_optimizer(o3c1e); call cint_del_optimizer(ogr)
      deallocate(fatm,fbas,fenv,fshls)
   end do

   do m = 1, 6
      print '(A,A12,A,I10,A,I0)', "  ", NAMES(m), "  values=", n(m), "  differing=", nb(m)
   end do
   if (sum(nb) > 0) stop 1
   print '(A)', "  RESULT: PASS (optimizer changes nothing)"

contains

   subroutine both(which, opt, nq_in)
      integer, intent(in) :: which, nq_in
      type(cint_opt_t), target, intent(in) :: opt
      integer :: mm
      b0(1:nq_in) = 0; b1(1:nq_in) = 0
      ws%opt => null()
      hv = run(which, b0)
      ws%opt => opt
      hv = run(which, b1)
      ws%opt => null()
      do mm = 1, nq_in
         n(which) = n(which) + 1
         if (b0(mm) /= b1(mm)) then
            nb(which) = nb(which) + 1
            if (nb(which) == 1) print '(A,A,A,I0,A,4I4,A,I0)', &
               "  FIRST ", trim(NAMES(which)), " basis=", ib, " shls", &
               fshls(0), fshls(1), fshls(2), fshls(3), " m=", mm
         end if
      end do
   end subroutine both

   logical function run(which, buf)
      integer, intent(in) :: which
      real(dp), intent(inout) :: buf(0:)
      select case (which)
      case (1); run = int1e_ovlp_sph (buf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
      case (2); run = int2e_sph      (buf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
      case (3); run = int3c2e_sph    (buf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
      case (4); run = int2c2e_sph    (buf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
      case (5); run = int3c1e_sph    (buf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
      case (6); run = int1e_grids_sph(buf, dims, fshls, fatm, 8, fbas, int(nbas), fenv, ws)
      end select
   end function run

end program opt_check
