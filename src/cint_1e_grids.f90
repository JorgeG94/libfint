!
! The grids one-electron integrals: <i| 1/|r-r_g| |j> evaluated at a list of
! points rather than at nuclei.
!
! Ported from src/g1e_grids.c and src/cint1e_grids.c.  This is what an
! electrostatic-potential calculation asks for, and it is the one family in
! the catalogue whose g array has a shape of its own.
!
! THE GRID INDEX IS INNERMOST AND BLOCKED.  Everywhere else in libcint the g
! array is (roots, i, j); here it is (grid, roots, i, j), with the grid index
! chunked into blocks of GRID_BLKSIZE so one block's worth of intermediates
! fits in cache.  That single change ripples through everything: the strides
! carry a factor of GRID_BLKSIZE, the Rys roots are solved once per grid point
! rather than once per shell pair, and the output has three indices where the
! other 1e integrals have two.  It is why this file exists instead of another
! arm on cint_1e.f90.
!
! Blocking also means the loop over grid blocks is outside the primitive
! loops, not inside: the whole shell pair is evaluated for one block of grids
! before moving to the next, so the pair data is computed once and reused.
!
module cint_1e_grids
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws, ws_ensure, ws_alloc_d, ws_alloc_i, &
                             ws_ensure_pd, ws_mark, ws_rewind, &
                             ws_opt_log_maxc, ws_opt_non0, ws_opt_idx
   use cint_screen,    only: pair_data, cint_set_pairdata, &
                             cint_log_max_pgto_coeff, cint_non0coeff_byshell
   use cint_g1e,       only: cint_init_int1e_envvars, cint_g1e_index_xyz, &
                             cint_common_fac_sp, cint_prim_to_ctr_0, cint_prim_to_ctr_1
   use cint_rys_roots, only: cint_rys_roots_lr, cint_rys_roots_sr
   use cint_bas,       only: NPRIM_OF, PTR_EXP, PTR_COEFF, ANG_OF, NCTR_OF
   use cint_cart2sph,  only: cint_c2s_ket_sph, RESULT_IN_GCART, RESULT_IN_GSPH
   use cint_opt,       only: cint_del_optimizer, opt_set_log_maxc, &
                             opt_set_non0coeff, opt_setij, opt_gen_idx, opt_finish
   implicit none
   private

   public :: cint_init_int1e_grids_envvars, cint_g0_1e_grids, cint_gout1e_grids
   public :: cint_1e_grids_drv, int1e_grids_cart, int1e_grids_sph
   ! cint_1e_grids_spinor sits on the same primitive loop; only the transform
   ! on the way out differs.
   public :: cint_1e_grids_loop
   public :: cint_nabla1i_grids, cint_nabla1j_grids, cint_x1i_grids, cint_x1j_grids
   public :: int1e_grids_optimizer
   public :: cint_all_1e_grids_optimizer

   ! src/cint_config.h.in.  Not a tuning knob to be re-picked here: the C's
   ! strides are built from it, and a different value would make the two
   ! layouts disagree while both remaining self-consistent.
   integer, parameter, public :: GRID_BLKSIZE = 104

   real(dp), parameter :: PI = 3.1415926535897932384626433832795029_dp
   real(dp), parameter :: EXPCUTOFF_SR = 40.0_dp

contains

   subroutine cint_init_int1e_grids_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      type(cint_env_vars), intent(inout) :: envs
      integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      integer :: dli, dlj, rys_order, nroots
      logical :: ibase
      real(dp) :: omega

      call cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)

      ! shls(2) and shls(3) are not shells here: they bracket the grid range,
      ! which is the C's way of passing two more integers through a fixed
      ! signature.
      envs%nfl = shls(3) - shls(2)                  ! ngrids
      envs%grids_off = int(env(PTR_GRIDS)) + shls(2) * 3

      envs%common_factor = 2.0_dp * PI &
         * cint_common_fac_sp(envs%i_l) * cint_common_fac_sp(envs%j_l)

      rys_order = envs%nrys_roots
      nroots = rys_order
      omega = env(PTR_RANGE_OMEGA)
      if (omega < 0.0_dp .and. rys_order <= 3) nroots = nroots * 2
      envs%rys_order = rys_order
      envs%nrys_roots = nroots

      ibase = envs%li_ceil > envs%lj_ceil
      if (ibase) then
         dli = envs%li_ceil + envs%lj_ceil + 1
         dlj = envs%lj_ceil + 1
         envs%rirj = envs%ri - envs%rj
      else
         dli = envs%li_ceil + 1
         dlj = envs%li_ceil + envs%lj_ceil + 1
         envs%rirj = envs%rj - envs%ri
      end if
      envs%g_stride_i = GRID_BLKSIZE * nroots
      envs%g_stride_j = GRID_BLKSIZE * nroots * dli
      envs%g_size     = GRID_BLKSIZE * nroots * dli * dlj
      envs%g_stride_k = envs%g_size
      envs%g_stride_l = envs%g_size
   end subroutine cint_init_int1e_grids_envvars

   ! The G tensor for one block of grid points.  gridsT is the block's
   ! coordinates, transposed so each Cartesian component is contiguous.
   function cint_g0_1e_grids(buf, og, ogrt, ou, orr, ot2, cutoff, envs) result(ok)
      real(dp), intent(inout) :: buf(0:)
      integer,  intent(in)    :: og, ogrt, ou, orr, ot2
      real(dp), intent(in)    :: cutoff
      type(cint_env_vars), intent(in) :: envs
      integer :: ok

      integer  :: ngrids, bgrids, nroots, rorder, gs, gx, gy, gz, w
      integer  :: n, i, j, ig, nmax, lj, di, dj, err
      integer  :: p0x, p0y, p0z, p1x, p1y, p1z, p2x, p2y, p2z
      integer  :: orgx, orgy, orgz
      real(dp) :: aij, aij2, x, fac1, omega, zeta, omega2, theta, sqrt_theta
      real(dp) :: a0, tau2, tau_theta, fac_theta, temp_cutoff
      real(dp) :: rijrx(0:2), rx(0:2)
      real(dp) :: ubuf(0:MXRYSROOTS-1), wbuf(0:MXRYSROOTS-1)

      ok = 0
      ngrids = envs%nfl
      bgrids = min(ngrids - envs%nfk, GRID_BLKSIZE)
      nroots = envs%nrys_roots
      gs = envs%g_size
      gx = og; gy = og + gs; gz = og + gs*2
      w = gz
      aij = envs%ai + envs%aj

      do i = 0, nroots - 1
         do ig = 0, bgrids - 1
            buf(gx + ig + GRID_BLKSIZE*i) = 1.0_dp
            buf(gy + ig + GRID_BLKSIZE*i) = 1.0_dp
         end do
      end do
      do ig = 0, bgrids - 1
         buf(orr + ig + GRID_BLKSIZE*0) = buf(ogrt + ig + GRID_BLKSIZE*0) - envs%rij(0)
         buf(orr + ig + GRID_BLKSIZE*1) = buf(ogrt + ig + GRID_BLKSIZE*1) - envs%rij(1)
         buf(orr + ig + GRID_BLKSIZE*2) = buf(ogrt + ig + GRID_BLKSIZE*2) - envs%rij(2)
      end do

      omega = envs%env(PTR_RANGE_OMEGA)
      zeta  = envs%env(PTR_RINV_ZETA)

      if (omega == 0.0_dp .and. zeta == 0.0_dp) then
         fac1 = envs%fac / aij
         do ig = 0, bgrids - 1
            x = aij * rgsq(buf, orr, ig)
            err = cint_rys_roots_lr(nroots, x, ubuf, wbuf)
            do i = 0, nroots - 1
               ! u stores t^2
               buf(ou + ig + GRID_BLKSIZE*i) = ubuf(i) / (ubuf(i) + 1.0_dp)
               buf(w  + ig + GRID_BLKSIZE*i) = wbuf(i) * fac1
            end do
         end do

      else if (omega < 0.0_dp) then
         ! short-range part of the range-separated Coulomb operator
         a0 = aij
         fac1 = envs%fac / aij
         if (zeta == 0.0_dp) then
            tau2 = 1.0_dp
            omega2 = omega * omega
            theta = omega2 / (omega2 + aij)
         else
            tau2 = zeta / (zeta + aij)
            a0 = a0 * tau2
            fac1 = fac1 * sqrt(tau2)
            omega2 = omega * omega
            theta = omega2 / (omega2 + a0)
         end if
         sqrt_theta = sqrt(theta)
         ! a very small erfc() gives ~0 weights, which upset the short-range
         ! root solver; the C caps the cutoff for that reason and so does this
         temp_cutoff = min(cutoff, EXPCUTOFF_SR)
         rorder = envs%rys_order
         do ig = 0, bgrids - 1
            x = a0 * rgsq(buf, orr, ig)
            if (theta * x > temp_cutoff) then
               do i = 0, nroots - 1
                  buf(ou + ig + GRID_BLKSIZE*i) = 0.0_dp
                  buf(w  + ig + GRID_BLKSIZE*i) = 0.0_dp
               end do
            else if (rorder == nroots) then
               err = cint_rys_roots_sr(nroots, x, sqrt_theta, ubuf, wbuf)
               do i = 0, nroots - 1
                  buf(ou + ig + GRID_BLKSIZE*i) = ubuf(i) / (ubuf(i) + 1.0_dp) * tau2
                  buf(w  + ig + GRID_BLKSIZE*i) = wbuf(i) * fac1
               end do
            else
               tau_theta = tau2 * theta
               fac_theta = fac1 * (-sqrt_theta)
               err = cint_rys_roots_lr(rorder, x, ubuf, wbuf)
               err = cint_rys_roots_lr(rorder, theta*x, ubuf(rorder:), wbuf(rorder:))
               do i = 0, rorder - 1
                  buf(ou + ig + GRID_BLKSIZE*i) = ubuf(i) / (ubuf(i) + 1.0_dp) * tau2
                  buf(w  + ig + GRID_BLKSIZE*i) = wbuf(i) * fac1
                  buf(ou + ig + GRID_BLKSIZE*(i+rorder)) = &
                     ubuf(i+rorder) / (ubuf(i+rorder) + 1.0_dp) * tau_theta
                  buf(w  + ig + GRID_BLKSIZE*(i+rorder)) = wbuf(i+rorder) * fac_theta
               end do
            end if
         end do

      else
         ! long-range part, or a Gaussian charge model
         a0 = aij
         fac1 = envs%fac / aij
         if (zeta == 0.0_dp) then
            omega2 = omega * omega
            theta = omega2 / (omega2 + aij)
         else if (omega == 0.0_dp) then
            theta = zeta / (zeta + aij)
         else
            omega2 = omega * omega
            theta = omega2*zeta / (omega2*zeta + (zeta + omega2)*aij)
         end if
         a0 = a0 * theta
         fac1 = fac1 * sqrt(theta)
         do ig = 0, bgrids - 1
            x = a0 * rgsq(buf, orr, ig)
            err = cint_rys_roots_lr(nroots, x, ubuf, wbuf)
            do i = 0, nroots - 1
               buf(ou + ig + GRID_BLKSIZE*i) = ubuf(i) / (ubuf(i) + 1.0_dp) * theta
               buf(w  + ig + GRID_BLKSIZE*i) = wbuf(i) * fac1
            end do
         end do
      end if

      nmax = envs%li_ceil + envs%lj_ceil
      if (nmax == 0) then
         ok = 1
         return
      end if

      if (envs%li_ceil > envs%lj_ceil) then
         lj = envs%lj_ceil
         di = envs%g_stride_i
         dj = envs%g_stride_j
         rx = envs%ri
      else
         lj = envs%li_ceil
         di = envs%g_stride_j
         dj = envs%g_stride_i
         rx = envs%rj
      end if
      rijrx = envs%rij - rx

      orgx = ot2 + GRID_BLKSIZE
      orgy = orgx + GRID_BLKSIZE
      orgz = orgy + GRID_BLKSIZE
      aij2 = 0.5_dp / aij

      do n = 0, nroots - 1
         p0x = gx + GRID_BLKSIZE*n; p0y = gy + GRID_BLKSIZE*n; p0z = gz + GRID_BLKSIZE*n
         p1x = p0x + di;            p1y = p0y + di;            p1z = p0z + di
         do ig = 0, bgrids - 1
            buf(orgx+ig) = rijrx(0) + buf(ou+ig+GRID_BLKSIZE*n) * buf(orr+ig+GRID_BLKSIZE*0)
            buf(orgy+ig) = rijrx(1) + buf(ou+ig+GRID_BLKSIZE*n) * buf(orr+ig+GRID_BLKSIZE*1)
            buf(orgz+ig) = rijrx(2) + buf(ou+ig+GRID_BLKSIZE*n) * buf(orr+ig+GRID_BLKSIZE*2)
            buf(p1x+ig) = buf(orgx+ig) * buf(p0x+ig)
            buf(p1y+ig) = buf(orgy+ig) * buf(p0y+ig)
            buf(p1z+ig) = buf(orgz+ig) * buf(p0z+ig)
         end do
         do ig = 0, bgrids - 1
            buf(ot2+ig) = aij2 * (1.0_dp - buf(ou+ig+GRID_BLKSIZE*n))
         end do
         do i = 1, nmax - 1
            p0x = gx + GRID_BLKSIZE*n + i*di
            p0y = gy + GRID_BLKSIZE*n + i*di
            p0z = gz + GRID_BLKSIZE*n + i*di
            p1x = p0x + di; p1y = p0y + di; p1z = p0z + di
            p2x = p0x - di; p2y = p0y - di; p2z = p0z - di
            do ig = 0, bgrids - 1
               buf(p1x+ig) = i * buf(ot2+ig) * buf(p2x+ig) + buf(orgx+ig) * buf(p0x+ig)
               buf(p1y+ig) = i * buf(ot2+ig) * buf(p2y+ig) + buf(orgy+ig) * buf(p0y+ig)
               buf(p1z+ig) = i * buf(ot2+ig) * buf(p2z+ig) + buf(orgz+ig) * buf(p0z+ig)
            end do
         end do
      end do

      do j = 1, lj
      do i = 0, nmax - j
         p0x = gx + j*dj + i*di; p0y = gy + j*dj + i*di; p0z = gz + j*dj + i*di
         p1x = p0x - dj;         p1y = p0y - dj;         p1z = p0z - dj
         p2x = p1x + di;         p2y = p1y + di;         p2z = p1z + di
         do n = 0, nroots - 1
         do ig = 0, bgrids - 1
            buf(p0x+ig+GRID_BLKSIZE*n) = buf(p2x+ig+GRID_BLKSIZE*n) &
                                       + envs%rirj(0) * buf(p1x+ig+GRID_BLKSIZE*n)
            buf(p0y+ig+GRID_BLKSIZE*n) = buf(p2y+ig+GRID_BLKSIZE*n) &
                                       + envs%rirj(1) * buf(p1y+ig+GRID_BLKSIZE*n)
            buf(p0z+ig+GRID_BLKSIZE*n) = buf(p2z+ig+GRID_BLKSIZE*n) &
                                       + envs%rirj(2) * buf(p1z+ig+GRID_BLKSIZE*n)
         end do
         end do
      end do
      end do
      ok = 1
   end function cint_g0_1e_grids

   pure function rgsq(buf, o, ig) result(r)
      real(dp), intent(in) :: buf(0:)
      integer,  intent(in) :: o, ig
      real(dp) :: r
      r = buf(o+ig+GRID_BLKSIZE*0)**2 + buf(o+ig+GRID_BLKSIZE*1)**2 &
        + buf(o+ig+GRID_BLKSIZE*2)**2
   end function rgsq

   ! Contract the three Cartesian factors over the Rys roots, for every grid
   ! point in the block.
   subroutine cint_gout1e_grids(gout, g, idx, envs, gout_empty)
      real(dp), intent(inout) :: gout(0:*)
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: idx(0:*)
      type(cint_env_vars), intent(in) :: envs
      integer,  intent(in)    :: gout_empty
      integer  :: ngrids, bgrids, nroots, nf, n, i, ig, ix, iy, iz
      real(dp) :: s(0:GRID_BLKSIZE-1)

      ngrids = envs%nfl
      bgrids = min(ngrids - envs%nfk, GRID_BLKSIZE)
      nroots = envs%nrys_roots
      nf = envs%nf

      do n = 0, nf - 1
         ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
         s(0:bgrids-1) = 0.0_dp
         do i = 0, nroots - 1
            do ig = 0, bgrids - 1
               s(ig) = s(ig) + g(ix+ig+GRID_BLKSIZE*i) * g(iy+ig+GRID_BLKSIZE*i) &
                                                       * g(iz+ig+GRID_BLKSIZE*i)
            end do
         end do
         if (gout_empty /= 0) then
            do ig = 0, bgrids - 1
               gout(ig + bgrids*n) = s(ig)
            end do
         else
            do ig = 0, bgrids - 1
               gout(ig + bgrids*n) = gout(ig + bgrids*n) + s(ig)
            end do
         end if
      end do
   end subroutine cint_gout1e_grids

   ! ---- the loop -------------------------------------------------------

   function cint_1e_grids_loop(gctr, gctroff, envs, ws) result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer,  intent(in)    :: gctroff
      type(cint_env_vars), intent(inout) :: envs
      ! TARGET for the block pointers bound below
      type(cint_ws), intent(inout), target :: ws
      logical :: has_value

      ! Hoisted as in the other loops.  goutp is rebound inside the grid-block
      ! loop, not once: with a single component the contraction writes
      ! straight into the caller's array and ogout moves with the block.
      real(dp), pointer, contiguous :: aip(:), ajp(:), cip(:), cjp(:)
      type(pair_data), pointer, contiguous :: pdp(:)
      real(dp), pointer, contiguous :: gp(:), goutp(:)
      integer,  pointer, contiguous :: idxp(:), nonp(:)
      integer :: i_sh, j_sh, i_ctr, j_ctr, i_prim, j_prim, nf, n_comp, ngrids
      logical :: use_opt
      logical :: idx_cached
      integer :: ai_o, aj_o, ci_o, cj_o, nc, leng, lenj, leni, len0
      integer :: oidx, olm, onon, osrt, og, ogout, ogctri, ogctrj, ogrt
      integer :: ou, orr, ot2, ip, jp, i, k, gofs, bgrids
      real(dp) :: expcutoff, rr_ij, fac1i, fac1j, cutoff
      logical  :: all_empty
      ! gempty, iempty, jempty, aliased as in the C
      integer, parameter :: E_G = 0, E_I = 1, E_J = 2
      logical :: flag(0:2)
      integer :: ig_, ii, ij

      has_value = .false.
      i_sh = envs%shls(0); j_sh = envs%shls(1)
      i_ctr = envs%x_ctr(0); j_ctr = envs%x_ctr(1)
      i_prim = bas_of(envs, NPRIM_OF, i_sh)
      j_prim = bas_of(envs, NPRIM_OF, j_sh)
      ai_o = bas_of(envs, PTR_EXP, i_sh);   aj_o = bas_of(envs, PTR_EXP, j_sh)
      ci_o = bas_of(envs, PTR_COEFF, i_sh); cj_o = bas_of(envs, PTR_COEFF, j_sh)
      nf = envs%nf
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor
      ngrids = envs%nfl
      expcutoff = envs%expcutoff
      rr_ij = sum(envs%rirj * envs%rirj)

      use_opt = .false.
      if (associated(ws%opt)) use_opt = cint_opt_usable(ws%opt, 2, envs)

      call ws_alloc_d(ws, i_prim + j_prim, olm)
      if (use_opt) then
         call ws_opt_log_maxc(ws, olm,          i_sh, i_prim)
         call ws_opt_log_maxc(ws, olm + i_prim, j_sh, j_prim)
      else
      call cint_log_max_pgto_coeff(ws%d(olm:),        envs%env, ci_o, i_prim, i_ctr)
      call cint_log_max_pgto_coeff(ws%d(olm+i_prim:), envs%env, cj_o, j_prim, j_ctr)
      end if
      call ws_ensure_pd(ws, i_prim*j_prim)
      if (cint_set_pairdata(ws%p, envs%env, ai_o, envs%env, aj_o, envs%ri, envs%rj, &
                            ws%d(olm:), ws%d(olm+i_prim:), envs%li_ceil, envs%lj_ceil, &
                            i_prim, j_prim, rr_ij, expcutoff, &
                            envs%env(PTR_RANGE_OMEGA))) then
         return
      end if

      call ws_alloc_i(ws, nf*3, oidx)
      ! NOT `use_opt .and. ws_opt_idx(...)`: Fortran does not
      ! guarantee short-circuit evaluation, so the compiler is free to
      ! call ws_opt_idx even when use_opt is false -- and that
      ! dereferences ws%opt, which is null precisely when use_opt is
      ! false.  gfortran short-circuits at -O2 and does not at -O0, so
      ! this segfaulted only in a debug build.
      idx_cached = .false.
      if (use_opt) idx_cached = ws_opt_idx(ws, oidx, cint_opt_idx_key(envs, 2), nf*3)
      if (.not. idx_cached) then
         call cint_g1e_index_xyz(ws%i(oidx:), envs)
      end if
      call ws_alloc_i(ws, i_prim + j_prim, onon)
      call ws_alloc_i(ws, i_prim*i_ctr + j_prim*j_ctr, osrt)
      if (use_opt) then
         call ws_opt_non0(ws, onon, osrt, i_sh, i_prim, i_ctr, 0, 0)
         call ws_opt_non0(ws, onon, osrt, j_sh, j_prim, j_ctr, i_prim, i_prim*i_ctr)
      else
      call cint_non0coeff_byshell(ws%i(osrt:), ws%i(onon:), envs%env, ci_o, i_prim, i_ctr)
      call cint_non0coeff_byshell(ws%i(osrt+i_prim*i_ctr:), ws%i(onon+i_prim:), &
                                  envs%env, cj_o, j_prim, j_ctr)
      end if

      nc = i_ctr * j_ctr
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenj = GRID_BLKSIZE * nf * nc * n_comp
      leni = GRID_BLKSIZE * nf * i_ctr * n_comp
      len0 = GRID_BLKSIZE * nf * n_comp

      call ws_alloc_d(ws, GRID_BLKSIZE*3, ogrt)          ! the transposed block
      call ws_alloc_d(ws, GRID_BLKSIZE*envs%nrys_roots, ou)
      call ws_alloc_d(ws, GRID_BLKSIZE*3, orr)
      call ws_alloc_d(ws, GRID_BLKSIZE*4, ot2)
      call ws_alloc_d(ws, leng, og)

      flag = .true.
      ig_ = E_G; ii = E_I; ij = E_J
      if (n_comp == 1) then
         ogctrj = gctroff
      else
         call ws_alloc_d(ws, lenj, ogctrj)
      end if
      if (j_ctr == 1) then
         ogctri = ogctrj; ii = ij
      else
         call ws_alloc_d(ws, leni, ogctri)
      end if
      if (i_ctr == 1) then
         ogout = ogctri; ig_ = ii
      else
         call ws_alloc_d(ws, len0, ogout)
      end if

      aip(0:) => envs%env(ai_o:); ajp(0:) => envs%env(aj_o:)
      cip(0:) => envs%env(ci_o:); cjp(0:) => envs%env(cj_o:)
      pdp(0:)  => ws%p
      gp(0:)   => ws%d(og:og+leng-1)
      idxp(0:) => ws%i(oidx:oidx+nf*3-1)
      nonp(0:) => ws%i(onon:onon+i_prim+j_prim-1)

      all_empty = .true.
      do gofs = 0, ngrids - 1, GRID_BLKSIZE
         envs%nfk = gofs                      ! grids_offset
         bgrids = min(ngrids - gofs, GRID_BLKSIZE)
         do i = 0, bgrids - 1
            ws%d(ogrt + i + GRID_BLKSIZE*0) = envs%env(envs%grids_off + (gofs+i)*3 + 0)
            ws%d(ogrt + i + GRID_BLKSIZE*1) = envs%env(envs%grids_off + (gofs+i)*3 + 1)
            ws%d(ogrt + i + GRID_BLKSIZE*2) = envs%env(envs%grids_off + (gofs+i)*3 + 2)
         end do

         flag = .true.
         ! With one component the contraction writes straight into the
         ! caller's array, at this block's offset.
         if (n_comp == 1) then
            ogctrj = gctroff + gofs * nf * nc
            if (j_ctr == 1) ogctri = ogctrj
            if (i_ctr == 1) ogout  = ogctri
         end if
         goutp(0:) => ws%d(ogout:ogout+len0-1)

         do jp = 0, j_prim - 1
            envs%aj = ajp(jp)
            if (j_ctr == 1) then
               fac1j = envs%common_factor * cjp(jp)
            else
               fac1j = envs%common_factor
               flag(ii) = .true.
            end if
            do ip = 0, i_prim - 1
               k = jp * i_prim + ip
               if (pdp(k)%cceij > expcutoff) cycle
               envs%ai = aip(ip)
               cutoff = expcutoff - pdp(k)%cceij
               envs%rij = pdp(k)%rij
               if (i_ctr == 1) then
                  fac1i = fac1j * cip(ip) * pdp(k)%eij
               else
                  fac1i = fac1j * pdp(k)%eij
               end if
               envs%fac = fac1i
               i = cint_g0_1e_grids(ws%d, og, ogrt, ou, orr, ot2, cutoff, envs)
               ! Through the procedure pointer, not to cint_gout1e_grids
               ! directly.  Calling the default kernel here makes every
               ! grids integral compute the plain one -- which is the D6
               ! bug, in a new file, found the same way.
               call envs%f_gout(goutp, gp, idxp, envs, merge(1, 0, flag(ig_)))
               if (i_ctr > 1) then
                  if (flag(ii)) then
                     call cint_prim_to_ctr_0(ws%d, ogctri, ogout, envs%env, ci_o+ip, &
                                             bgrids*nf*n_comp, i_prim, i_ctr)
                  else
                     call cint_prim_to_ctr_1(ws%d, ogctri, ogout, envs%env, ci_o+ip, &
                                             bgrids*nf*n_comp, i_prim, i_ctr, &
                                             nonp(ip), ws%i, osrt + ip*i_ctr)
                  end if
               end if
               flag(ii) = .false.
            end do
            if (.not. flag(ii)) then
               if (j_ctr > 1) then
                  if (flag(ij)) then
                     call cint_prim_to_ctr_0(ws%d, ogctrj, ogctri, envs%env, cj_o+jp, &
                                             bgrids*nf*i_ctr*n_comp, j_prim, j_ctr)
                  else
                     call cint_prim_to_ctr_1(ws%d, ogctrj, ogctri, envs%env, cj_o+jp, &
                                             bgrids*nf*i_ctr*n_comp, j_prim, j_ctr, &
                                             nonp(i_prim+jp), ws%i, &
                                             osrt + i_prim*i_ctr + jp*j_ctr)
                  end if
               end if
               flag(ij) = .false.
            end if
         end do

         if (n_comp > 1 .and. .not. flag(ij)) then
            call transpose_comps(gctr, gctroff + gofs*nf*nc, ws%d, ogctrj, &
                                 bgrids, nf*nc, ngrids, n_comp)
         end if
         all_empty = all_empty .and. flag(ij)
      end do
      has_value = .not. all_empty
   end function cint_1e_grids_loop


   ! The contraction leaves the components interleaved per (i,j); the caller
   ! wants them as separate planes.
   pure subroutine transpose_comps(gctr, goff, src, soff, bgrids, dij, ngrids, n_comp)
      real(dp), intent(inout) :: gctr(0:)
      real(dp), intent(in)    :: src(0:)
      integer,  intent(in)    :: goff, soff, bgrids, dij, ngrids, n_comp
      integer :: n, ic, ig, po, ps
      do ic = 0, n_comp - 1
         po = goff + ic * dij * ngrids
         do n = 0, dij - 1
            ps = soff + (n * n_comp + ic) * bgrids
            do ig = 0, bgrids - 1
               gctr(po + ig + n*bgrids) = src(ps + ig)
            end do
         end do
      end do
   end subroutine transpose_comps

   ! ---- output ---------------------------------------------------------

   pure subroutine dcopy_grids_ij(out, ooff, src, soff, ngrids, ni, mgrids, mi, mj)
      real(dp), intent(inout) :: out(0:)
      real(dp), intent(in)    :: src(0:)
      integer,  intent(in)    :: ooff, soff, ngrids, ni, mgrids, mi, mj
      integer :: i, j, m, ngi, mgi, ob, sb
      ngi = ngrids * ni
      mgi = mgrids * mi
      ob = ooff; sb = soff
      do j = 0, mj - 1
         do i = 0, mi - 1
            do m = 0, mgrids - 1
               out(ob + i*ngrids + m) = src(sb + i*mgrids + m)
            end do
         end do
         ob = ob + ngi
         sb = sb + mgi
      end do
   end subroutine dcopy_grids_ij

   subroutine apply_c2s_cart_1e_grids(out, buf, goff, dims, envs)
      real(dp), intent(inout) :: out(0:)
      real(dp), intent(in)    :: buf(0:)
      integer,  intent(in)    :: goff, dims(0:)
      type(cint_env_vars), intent(in) :: envs
      integer :: ngrids, i_ctr, j_ctr, ni, ng, nfi, nfj, nf, ofj
      integer :: ic, jc, gofs, bgrids, gb
      ngrids = envs%nfl
      i_ctr = envs%x_ctr(0); j_ctr = envs%x_ctr(1)
      ni = dims(0); ng = dims(2)
      nfi = envs%nfi; nfj = envs%nfj; nf = envs%nf
      ofj = ni * nfj
      gb = goff
      do gofs = 0, ngrids - 1, GRID_BLKSIZE
         bgrids = min(ngrids - gofs, GRID_BLKSIZE)
         do jc = 0, j_ctr - 1
         do ic = 0, i_ctr - 1
            call dcopy_grids_ij(out, ng*(ofj*jc + nfi*ic) + gofs, buf, gb, &
                                ng, ni, bgrids, nfi, nfj)
            gb = gb + bgrids * nf
         end do
         end do
      end do
   end subroutine apply_c2s_cart_1e_grids

   subroutine apply_c2s_sph_1e_grids(out, buf, goff, o1, o2, dims, envs)
      real(dp), intent(inout) :: out(0:)
      real(dp), intent(inout) :: buf(0:)
      integer,  intent(in)    :: goff, o1, o2, dims(0:)
      type(cint_env_vars), intent(in) :: envs
      integer :: ngrids, i_l, j_l, i_ctr, j_ctr, di, dj, ni, ng, ofj, nfi, nf
      integer :: ic, jc, gofs, bgrids, bgdi, bgnfi, gb, loc, cur, n, r, buflen
      ngrids = envs%nfl
      i_l = envs%i_l; j_l = envs%j_l
      i_ctr = envs%x_ctr(0); j_ctr = envs%x_ctr(1)
      di = i_l*2+1; dj = j_l*2+1
      ni = dims(0); ng = dims(2)
      ofj = ni * dj
      nfi = envs%nfi; nf = envs%nf
      buflen = GRID_BLKSIZE * nfi * dj

      gb = goff
      do gofs = 0, ngrids - 1, GRID_BLKSIZE
         bgrids = min(ngrids - gofs, GRID_BLKSIZE)
         bgdi = bgrids * di
         bgnfi = bgrids * nfi
         do jc = 0, j_ctr - 1
         do ic = 0, i_ctr - 1
            ! the ket first, then the bra -- and each is the identity for
            ! l <= 1, which is what the loc sentinel reports
            loc = cint_c2s_ket_sph(buf(o1:), buf(gb:), bgnfi, bgnfi, j_l)
            if (loc == RESULT_IN_GCART) then
               buf(o1:o1+buflen-1) = buf(gb:gb+buflen-1)
            end if
            if (i_l <= 1) then
               cur = o1
            else
               do n = 0, dj - 1
                  r = cint_c2s_ket_sph(buf(o2 + n*bgdi:), buf(o1 + n*bgnfi:), &
                                       bgrids, bgrids, i_l)
               end do
               cur = o2
            end if
            call dcopy_grids_ij(out, ng*(ofj*jc + di*ic) + gofs, buf, cur, &
                                ng, ni, bgrids, di, dj)
            gb = gb + bgrids * nf
         end do
         end do
      end do
   end subroutine apply_c2s_sph_1e_grids

   pure subroutine c2s_grids_dset0(out, dims, counts)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), counts(0:)
      integer :: i, j, m, ni, ng
      ni = dims(0); ng = dims(2)
      do j = 0, counts(1) - 1
      do i = 0, counts(0) - 1
      do m = 0, counts(2) - 1
         out(ng*(j*ni + i) + m) = 0.0_dp
      end do
      end do
      end do
   end subroutine c2s_grids_dset0

   ! ---- driver ---------------------------------------------------------

   function cint_1e_grids_drv(out, dims, envs, ws, c2s_kind) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      integer,  intent(in)    :: c2s_kind
      logical :: has_value

      integer :: ngrids, nf, nc, n_comp, ogctr, o1, o2, n, nout, counts(0:3)
      integer :: nd, ni, ip, jp, leng, len0, leni, lenj, buflen, dmark, imark

      ngrids = envs%nfl
      nf = envs%nf
      nc = ngrids * nf * envs%x_ctr(0) * envs%x_ctr(1)
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor

      ip = bas_of(envs, NPRIM_OF, envs%shls(0))
      jp = bas_of(envs, NPRIM_OF, envs%shls(1))
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      len0 = GRID_BLKSIZE * nf * n_comp
      leni = len0 * envs%x_ctr(0)
      lenj = leni * envs%x_ctr(1)
      buflen = GRID_BLKSIZE * envs%nfi * (2*envs%j_l + 1)
      nd = nc*n_comp + ip + jp + leng + len0 + leni + lenj &
         + GRID_BLKSIZE*(10 + envs%nrys_roots) + 2*buflen + 128
      ni = nf*3 + ip + jp + ip*envs%x_ctr(0) + jp*envs%x_ctr(1) + 64
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*n_comp, ogctr)
      ! NOT zeroed.  No C driver zeroes gctr -- CINT1e_drv, CINT3c2e_drv,
      ! CINT2c2e_drv, CINT3c1e_drv and CINT1e_grids_drv all just take the
      ! buffer and thread an empty flag through the loop, so the first
      ! contraction sets rather than accumulates.  The port threads the same
      ! flag, so the memset was pure extra writes.

      has_value = cint_1e_grids_loop(ws%d, ogctr, envs, ws)

      select case (c2s_kind)
      case (C2S_SPH_1E)
         counts(0) = (envs%i_l*2+1) * envs%x_ctr(0)
         counts(1) = (envs%j_l*2+1) * envs%x_ctr(1)
      case default
         counts(0) = envs%nfi * envs%x_ctr(0)
         counts(1) = envs%nfj * envs%x_ctr(1)
      end select
      counts(2) = ngrids
      counts(3) = 1
      nout = dims(0) * dims(1) * dims(2)

      if (has_value) then
         call ws_mark(ws, dmark, imark)
         do n = 0, n_comp - 1
            call ws_rewind(ws, dmark, imark)
            select case (c2s_kind)
            case (C2S_SPH_1E)
               call ws_alloc_d(ws, buflen, o1)
               call ws_alloc_d(ws, buflen, o2)
               call apply_c2s_sph_1e_grids(out(nout*n:), ws%d, ogctr + nc*n, o1, o2, &
                                           dims, envs)
            case default
               call apply_c2s_cart_1e_grids(out(nout*n:), ws%d, ogctr + nc*n, dims, envs)
            end select
         end do
      else
         do n = 0, n_comp - 1
            call c2s_grids_dset0(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_1e_grids_drv

   ! ---- entry points ---------------------------------------------------

   function int1e_grids_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int1e_grids_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout1e_grids
      hv = cint_1e_grids_drv(out, dims, envs, ws, C2S_CART_1E)
   end function int1e_grids_cart

   function int1e_grids_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int1e_grids_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout1e_grids
      hv = cint_1e_grids_drv(out, dims, envs, ws, C2S_SPH_1E)
   end function int1e_grids_sph


   ! ---- the grids recursion operations ---------------------------------
   !
   ! What the G1E_GRIDS_* macros expand to.  Same shape as their plain 1e
   ! counterparts in cint_g1e.f90, with the grid index innermost: every
   ! innermost loop runs over a block of grid points rather than over Rys
   ! roots, and `ptr` carries the n*GRID_BLKSIZE the roots contribute.
   !
   ! One array, two offsets, for the aliasing reason the rest of this port
   ! carries.

   subroutine cint_nabla1i_grids(g, foff, goff, li, lj, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj
      type(cint_env_vars), intent(in) :: envs
      integer  :: ngrids, bgrids, nroots, di, dj, gs, i, j, n, ig, ptr
      integer  :: fx, fy, fz, gx, gy, gz
      real(dp) :: ai2
      ngrids = envs%nfl
      bgrids = min(ngrids - envs%nfk, GRID_BLKSIZE)
      nroots = envs%nrys_roots
      di = envs%g_stride_i; dj = envs%g_stride_j; gs = envs%g_size
      ai2 = -2.0_dp * envs%ai
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2

      do j = 0, lj
         do n = 0, nroots - 1
            ptr = dj*j + n*GRID_BLKSIZE
            do ig = ptr, ptr + bgrids - 1
               g(fx+ig) = ai2 * g(gx+ig+di)
               g(fy+ig) = ai2 * g(gy+ig+di)
               g(fz+ig) = ai2 * g(gz+ig+di)
            end do
         end do
         do i = 1, li
         do n = 0, nroots - 1
            ptr = dj*j + di*i + n*GRID_BLKSIZE
            do ig = ptr, ptr + bgrids - 1
               g(fx+ig) = i * g(gx+ig-di) + ai2 * g(gx+ig+di)
               g(fy+ig) = i * g(gy+ig-di) + ai2 * g(gy+ig+di)
               g(fz+ig) = i * g(gz+ig-di) + ai2 * g(gz+ig+di)
            end do
         end do
         end do
      end do
   end subroutine cint_nabla1i_grids

   subroutine cint_nabla1j_grids(g, foff, goff, li, lj, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj
      type(cint_env_vars), intent(in) :: envs
      integer  :: ngrids, bgrids, nroots, di, dj, gs, i, j, n, ig, ptr
      integer  :: fx, fy, fz, gx, gy, gz
      real(dp) :: aj2
      ngrids = envs%nfl
      bgrids = min(ngrids - envs%nfk, GRID_BLKSIZE)
      nroots = envs%nrys_roots
      di = envs%g_stride_i; dj = envs%g_stride_j; gs = envs%g_size
      aj2 = -2.0_dp * envs%aj
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2

      do i = 0, li
      do n = 0, nroots - 1
         ptr = di*i + n*GRID_BLKSIZE
         do ig = ptr, ptr + bgrids - 1
            g(fx+ig) = aj2 * g(gx+ig+dj)
            g(fy+ig) = aj2 * g(gy+ig+dj)
            g(fz+ig) = aj2 * g(gz+ig+dj)
         end do
      end do
      end do
      do j = 1, lj
      do i = 0, li
      do n = 0, nroots - 1
         ptr = dj*j + di*i + n*GRID_BLKSIZE
         do ig = ptr, ptr + bgrids - 1
            g(fx+ig) = j * g(gx+ig-dj) + aj2 * g(gx+ig+dj)
            g(fy+ig) = j * g(gy+ig-dj) + aj2 * g(gy+ig+dj)
            g(fz+ig) = j * g(gz+ig-dj) + aj2 * g(gz+ig+dj)
         end do
      end do
      end do
      end do
   end subroutine cint_nabla1j_grids

   subroutine cint_x1i_grids(g, foff, goff, r, li, lj, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj
      real(dp), intent(in)    :: r(0:2)
      type(cint_env_vars), intent(in) :: envs
      call x1_grids(g, foff, goff, r, li, lj, envs, envs%g_stride_i)
   end subroutine cint_x1i_grids

   subroutine cint_x1j_grids(g, foff, goff, r, li, lj, envs)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj
      real(dp), intent(in)    :: r(0:2)
      type(cint_env_vars), intent(in) :: envs
      call x1_grids(g, foff, goff, r, li, lj, envs, envs%g_stride_j)
   end subroutine cint_x1j_grids

   subroutine x1_grids(g, foff, goff, r, li, lj, envs, stride)
      real(dp), intent(inout) :: g(0:)
      integer,  intent(in)    :: foff, goff, li, lj, stride
      real(dp), intent(in)    :: r(0:2)
      type(cint_env_vars), intent(in) :: envs
      integer :: ngrids, bgrids, nroots, di, dj, gs, i, j, n, ig, ptr
      integer :: fx, fy, fz, gx, gy, gz
      ngrids = envs%nfl
      bgrids = min(ngrids - envs%nfk, GRID_BLKSIZE)
      nroots = envs%nrys_roots
      di = envs%g_stride_i; dj = envs%g_stride_j; gs = envs%g_size
      fx = foff; fy = foff + gs; fz = foff + gs*2
      gx = goff; gy = goff + gs; gz = goff + gs*2

      do j = 0, lj
      do i = 0, li
      do n = 0, nroots - 1
         ptr = dj*j + di*i + n*GRID_BLKSIZE
         do ig = ptr, ptr + bgrids - 1
            g(fx+ig) = g(gx+ig+stride) + r(0) * g(gx+ig)
            g(fy+ig) = g(gy+ig+stride) + r(1) * g(gy+ig)
            g(fz+ig) = g(gz+ig+stride) + r(2) * g(gz+ig)
         end do
      end do
      end do
      end do
   end subroutine x1_grids

   ! ---- optimizers ------------------------------------------------------

   subroutine int1e_grids_optimizer(opt, atm, natm, bas, nbas, env)
      ! The C's is `*opt = NULL`: this family has no optimizer upstream.
      ! Clearing rather than ignoring, so a stale optimizer left on the
      ! workspace cannot be picked up by the next call.
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
   end subroutine int1e_grids_optimizer


   ! ---- the optimizer builder for this arity ----------------------------
   !
   ! It lives here rather than in cint_opt because it names this module's
   ! envs init, and cint_opt has to stay underneath every arity.

   subroutine cint_all_1e_grids_optimizer(opt, ng, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: ng(0:), natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
      call opt_set_log_maxc(opt, bas, nbas, env)
      call opt_set_non0coeff(opt, bas, nbas, env)
      call opt_gen_idx(opt, cint_init_int1e_grids_envvars, cint_g1e_index_xyz, 2, ANG_MAX, &
                       ng, atm, natm, bas, nbas, env)
      call opt_finish(opt, 2, ng, nbas)
   end subroutine cint_all_1e_grids_optimizer


end module cint_1e_grids
