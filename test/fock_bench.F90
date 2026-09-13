!
! A Fock build: the four two-electron paths compared where it matters, and
! the only check that runs any of them under OpenMP.
!
! Schwarz-screened, threaded over shell quartets, with the same density and
! the same screening for every path, so the only variable is which routine
! evaluates the quartet.  It is a test as well as a benchmark -- the paths
! must agree on the assembled matrix, which is a different question from
! the per-quartet checks: it also says the drivers are re-entrant, keep no
! state between calls, and survive a threadprivate workspace.
!
! The per-quartet harnesses in this project mispredicted the Fock-build
! result twice, in both directions, which is why this exists.
!
! Water cluster, STO-3G (real data: an L shell on every oxygen, which is
! where the rotated-axis path is strongest), optionally with a d shell on
! the oxygens so that the d-touching classes appear at all.
!
!   ./fock_bench <nwater> <d: 0|1> <path: 0=all>
!
! The workspace is threadprivate, exactly as include/libfint.f90 declares
! its own: an OpenMP `private` copy of a derived type with allocatable
! components is undefined on entry, and the integral drivers bump-allocate
! from it.
module fb_ws
   use cint_workspace, only: cint_ws
   implicit none
   type(cint_ws), save :: ws
   !$omp threadprivate(ws)
end module fb_ws

program fock_bench
   use fb_ws, only: ws
   use iso_c_binding
   use cint_const, only: dp
   use cint_bas, only: ATM_SLOTS, BAS_SLOTS, ANG_OF, NCTR_OF, NPRIM_OF, KAPPA_OF, &
                       PTR_EXP, PTR_COEFF, ATOM_OF, PTR_COORD, CHARGE_OF, &
                       KAPPA_SP_SHELL, cint_gto_norm, cint_cgto_spheric, cint_bas_is_sp
   use cint_envs, only: PTR_ENV_START
   use cint_workspace, only: cint_ws
   use cint_2e, only: int2e_sph
   use cint_rotaxis_2e, only: int2e_rotaxis_sph, rotaxis_supported
   use cint_hgp_2e, only: int2e_hgp_sph, hgp_supported
   !$ use omp_lib
   implicit none

   integer, parameter :: MAXA = 400, MAXB = 2000, MAXE = 200000
   integer  :: atm(0:ATM_SLOTS*MAXA-1), bas(0:BAS_SLOTS*MAXB-1), natm, nbas, envoff
   real(dp) :: env(0:MAXE-1)
   integer  :: ao_loc(0:MAXB), nao
   real(dp), allocatable :: dm(:,:), fk(:,:,:), q(:), fref(:,:)
   integer,  allocatable :: pi(:), pj(:)
   integer   :: npair, nwat, dpol, only, ipath, i, j, nthr
   real(dp)  :: t0, t1, tsec(0:3), dmax, scal, rel, worst
   character(len=16) :: arg
   character(len=8), parameter :: pname(0:3) = ["rys     ", "rotaxis ", "hgp     ", "hybrid  "]

   nwat = 12; dpol = 1; only = 0; worst = 0.0_dp
   if (command_argument_count() >= 1) then
      call get_command_argument(1, arg); read(arg,*) nwat
   end if
   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg); read(arg,*) dpol
   end if
   if (command_argument_count() >= 3) then
      call get_command_argument(3, arg); read(arg,*) only
   end if

   call build(nwat, dpol)
   nthr = 1
   !$ nthr = omp_get_max_threads()
   print '(A,I0,A,I0,A,I0,A,I0,A,I0)', "  waters ", nwat, "  atoms ", natm, &
        "  shells ", nbas, "  basis functions ", nao, "  threads ", nthr

   allocate(dm(nao,nao), fref(nao,nao))
   call fake_density()
   dmax = maxval(abs(dm))
   call schwarz()
   print '(A,I0,A,I0,A)', "  shell pairs ", npair, " of ", nbas*(nbas+1)/2, " survive screening"

   do ipath = 0, 3
      if (only /= 0 .and. ipath /= only-1) cycle
      t0 = wall()
      call fock(ipath)
      t1 = wall()
      tsec(ipath) = t1 - t0
      if (ipath == 0) then
         fref = fk(:,:,0)
         scal = max(maxval(abs(fref)), 1.0e-30_dp)
         print '(A,A,F9.3,A)', "  ", pname(ipath), tsec(ipath), " s"
      else
         rel = maxval(abs(fk(:,:,0) - fref))/scal
         worst = max(worst, rel)
         print '(A,A,F9.3,A,F6.2,A,ES10.2)', "  ", pname(ipath), tsec(ipath), &
              " s   ", tsec(0)/tsec(ipath), "x rys   scaled |dF| ", rel
      end if
      deallocate(fk)
   end do

   ! Two exact algorithms round differently, so this is a tolerance and not
   ! bit identity -- the same argument rotaxis_check and hgp_check make,
   ! here on a matrix summed over millions of integrals rather than on one
   ! quartet.
   if (only == 0) then
      if (worst > 1.0e-11_dp) then
         print '(A,ES10.2)', "  RESULT: FAIL, worst scaled difference ", worst
         stop 1
      end if
      print '(A,ES10.2,A)', "  RESULT: PASS (four paths agree to ", worst, " scaled)"
   end if

contains

   function wall() result(t)
      real(dp) :: t
      integer(8) :: c, r
      call system_clock(c, r)
      t = real(c,dp)/real(r,dp)
   end function wall

   ! ---- the Fock build -------------------------------------------------

   subroutine fock(path)
      integer, intent(in) :: path
      integer  :: ip_, kp_, ish, jsh, ksh, lsh, di, dj, dk, dl
      integer  :: a, b, c, d, ii, jj, kk, ll, shls(0:3), dims(0:3), tid
      real(dp) :: v, s4
      real(dp), allocatable :: buf(:)
      logical  :: hv

      allocate(fk(nao,nao,0:nthr-1))
      fk = 0.0_dp
      !$omp parallel default(shared) private(ip_,kp_,ish,jsh,ksh,lsh,di,dj,dk,dl, &
      !$omp   a,b,c,d,ii,jj,kk,ll,shls,dims,v,s4,buf,hv,tid)
      tid = 0
      !$ tid = omp_get_thread_num()
      allocate(buf(0:15**4-1))
      !$omp do schedule(dynamic, 4)
      do ip_ = 1, npair
         do kp_ = 1, ip_
            if (q(ip_)*q(kp_)*dmax < 1.0e-10_dp) cycle
            ish = pi(ip_); jsh = pj(ip_); ksh = pi(kp_); lsh = pj(kp_)
            shls = [ish, jsh, ksh, lsh]
            di = ao_loc(ish+1)-ao_loc(ish); dj = ao_loc(jsh+1)-ao_loc(jsh)
            dk = ao_loc(ksh+1)-ao_loc(ksh); dl = ao_loc(lsh+1)-ao_loc(lsh)
            dims = [di, dj, dk, dl]
            hv = eval(path, buf, dims, shls)
            if (.not. hv) cycle
            ! the eight-fold symmetry factor for this quartet
            s4 = 1.0_dp
            if (ish /= jsh) s4 = s4*2.0_dp
            if (ksh /= lsh) s4 = s4*2.0_dp
            if (ip_ /= kp_) s4 = s4*2.0_dp
            do d = 0, dl-1
            do c = 0, dk-1
            do b = 0, dj-1
            do a = 0, di-1
               v = buf(a + di*(b + dj*(c + dk*d)))*s4
               if (abs(v) < 1.0e-14_dp) cycle
               ii = ao_loc(ish)+a; jj = ao_loc(jsh)+b
               kk = ao_loc(ksh)+c; ll = ao_loc(lsh)+d
               fk(ii+1,jj+1,tid) = fk(ii+1,jj+1,tid) + v*dm(kk+1,ll+1)
               fk(ii+1,kk+1,tid) = fk(ii+1,kk+1,tid) - 0.25_dp*v*dm(jj+1,ll+1)
            end do
            end do
            end do
            end do
         end do
      end do
      !$omp end do
      deallocate(buf)
      !$omp end parallel
      do i = 1, nthr-1
         fk(:,:,0) = fk(:,:,0) + fk(:,:,i)
      end do
   end subroutine fock

   logical function eval(path, buf, dims, shls) result(hv)
      integer,  intent(in) :: path, dims(0:), shls(0:)
      real(dp), intent(inout) :: buf(0:)
      select case (path)
      case (0); hv = int2e_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
      case (1)
         if (rotaxis_supported(shls, bas)) then
            hv = int2e_rotaxis_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         else
            hv = int2e_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         end if
      case (2)
         if (hgp_supported(shls, bas)) then
            hv = int2e_hgp_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         else
            hv = int2e_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         end if
      case default
         if (rotaxis_supported(shls, bas)) then
            hv = int2e_rotaxis_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         else if (hgp_supported(shls, bas)) then
            hv = int2e_hgp_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         else
            hv = int2e_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
         end if
      end select
   end function eval

   ! ---- screening -------------------------------------------------------

   subroutine schwarz()
      integer  :: ish, jsh, di, dj, shls(0:3), dims(0:3), m, n
      real(dp) :: buf(15**4), qq
      logical  :: hv
      allocate(q(nbas*(nbas+1)/2), pi(nbas*(nbas+1)/2), pj(nbas*(nbas+1)/2))
      npair = 0
      do ish = 0, nbas-1
         do jsh = 0, ish
            di = ao_loc(ish+1)-ao_loc(ish); dj = ao_loc(jsh+1)-ao_loc(jsh)
            shls = [ish, jsh, ish, jsh]; dims = [di, dj, di, dj]
            n = di*dj*di*dj
            buf(1:n) = 0.0_dp
            hv = int2e_sph(buf, dims, shls, atm, natm, bas, nbas, env, ws)
            qq = 0.0_dp
            do m = 1, n
               qq = max(qq, abs(buf(m)))
            end do
            qq = sqrt(qq)
            if (qq > 1.0e-12_dp) then
               npair = npair + 1
               q(npair) = qq; pi(npair) = ish; pj(npair) = jsh
            end if
         end do
      end do
   end subroutine schwarz

   subroutine fake_density()
      integer :: a, b
      integer(8) :: s
      s = 12345_8
      do a = 1, nao
         do b = 1, a
            s = mod(s*6364136223846793005_8 + 1442695040888963407_8, 2147483647_8)
            dm(a,b) = real(s,dp)/2147483647.0_dp - 0.5_dp
            dm(b,a) = dm(a,b)
         end do
      end do
   end subroutine fake_density

   ! ---- the system ------------------------------------------------------

   subroutine build(nw, withd)
      integer, intent(in) :: nw, withd
      real(dp) :: o3(3), oc(3), h3(3), hc(3), l3(3), lc(3,2), d1(1), dc(1,1)
      real(dp) :: r(3), ang
      integer  :: w, nx, ia
      atm = 0; bas = 0; env = 0.0_dp; natm = 0; nbas = 0; envoff = PTR_ENV_START
      ! STO-3G
      o3 = [130.7093200_dp, 23.8088610_dp, 6.4436083_dp]
      oc = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
      l3 = [5.0331513_dp, 1.1695961_dp, 0.3803890_dp]
      lc(:,1) = [-0.09996723_dp, 0.39951283_dp, 0.70011547_dp]
      lc(:,2) = [0.15591627_dp, 0.60768372_dp, 0.39195739_dp]
      h3 = [3.42525091_dp, 0.62391373_dp, 0.16885540_dp]
      hc = [0.15432897_dp, 0.53532814_dp, 0.44463454_dp]
      d1 = [0.8_dp]; dc = reshape([1.0_dp], [1,1])
      nx = ceiling(real(nw,dp)**(1.0_dp/3.0_dp))
      do w = 0, nw-1
         r = [real(mod(w, nx),dp), real(mod(w/nx, nx),dp), real(w/(nx*nx),dp)]*5.4_dp
         ang = 0.37_dp*w
         ia = natm
         call put_atom(8, r)
         call put_atom(1, r + [1.43_dp*cos(ang), 1.43_dp*sin(ang), 0.30_dp])
         call put_atom(1, r + [-0.45_dp, 1.36_dp*cos(ang), -1.20_dp])
         call put_shell(ia,   0, 3, 1, o3, reshape(oc,[3,1]), .false.)
         call put_shell(ia,   1, 3, 1, l3, lc, .true.)
         if (withd /= 0) call put_shell(ia, 2, 1, 1, d1, dc, .false.)
         call put_shell(ia+1, 0, 3, 1, h3, reshape(hc,[3,1]), .false.)
         call put_shell(ia+2, 0, 3, 1, h3, reshape(hc,[3,1]), .false.)
      end do
      ao_loc(0) = 0
      do i = 0, nbas-1
         ao_loc(i+1) = ao_loc(i) + cint_cgto_spheric(i, bas)
      end do
      nao = ao_loc(nbas)
   end subroutine build

   subroutine put_atom(z, r)
      integer,  intent(in) :: z
      real(dp), intent(in) :: r(3)
      atm(ATM_SLOTS*natm + CHARGE_OF) = z
      atm(ATM_SLOTS*natm + PTR_COORD) = envoff
      env(envoff:envoff+2) = r; envoff = envoff + 3
      natm = natm + 1
   end subroutine put_atom

   subroutine put_shell(ia, l, np, nc, e, c, is_sp)
      integer,  intent(in) :: ia, l, np, nc
      real(dp), intent(in) :: e(:), c(:,:)
      logical,  intent(in) :: is_sp
      integer :: ncol, m, ip_, eptr, cptr, lm, b
      ncol = nc; if (is_sp) ncol = 2*nc
      eptr = envoff; env(eptr:eptr+np-1) = e(1:np); envoff = envoff + np
      cptr = envoff
      do m = 1, ncol
         lm = l
         if (is_sp) then
            lm = 1; if (m <= nc) lm = 0
         end if
         do ip_ = 1, np
            env(cptr + (m-1)*np + ip_-1) = c(ip_,m)*cint_gto_norm(lm, e(ip_))
         end do
      end do
      envoff = envoff + ncol*np
      b = BAS_SLOTS*nbas
      bas(b+ATOM_OF) = ia; bas(b+ANG_OF) = l; bas(b+NPRIM_OF) = np
      bas(b+NCTR_OF) = nc; bas(b+KAPPA_OF) = 0
      if (is_sp) bas(b+KAPPA_OF) = KAPPA_SP_SHELL
      bas(b+PTR_EXP) = eptr; bas(b+PTR_COEFF) = cptr
      nbas = nbas + 1
   end subroutine put_shell

end program fock_bench
