!
! The two-electron G tensor.
!
! Ported from src/g2e.c.  Three layers:
!
!   * cint_init_int2e_envvars   -- strides, angular momenta, which transfer
!   * cint_g0_2e_2d             -- the Rys 2D recursion in (m,n)
!   * cint_g0_{lj,kj,il,ik}2d_4d -- transfer the 2D result into the 4D
!                                   (i,k,l,j) layout the gout kernels index
!   * cint_g0_2e                -- Rys roots, the b/c coefficients, dispatch
!
! Every stride and every index expression is the C's, unchanged
! (PORT_TO_FORTRAN.md 3.6).  The `g` array is one flat block of three
! Cartesian components, each g_size long, and gx/gy/gz are offsets into it
! rather than the C's three pointers.
!
module cint_g2e
   use cint_const,  only: dp
   use cint_envs
   use cint_bas,    only: cint_cart_comp, cint_cart_comp_sp, cint_bas_is_sp, &
                          NF_SP, ANG_OF, NCTR_OF, ATOM_OF, PTR_COORD
   use cint_g1e,    only: cint_common_fac_sp, cint_g1e_index_xyz
   use cint_rys_roots, only: cint_rys_roots_lr, cint_rys_roots_sr
   use cint_g2e_unrolled, only: cint_g0_2e_2d4d_unrolled, cint_srg0_2e_2d4d_unrolled
   use cint_opt,       only: cint_del_optimizer, opt_set_log_maxc, &
                             opt_set_non0coeff, opt_setij, opt_gen_idx, opt_finish
   implicit none
   private

   public :: cint_init_int2e_envvars, cint_init_int3c2e_envvars
   public :: cint_init_int2c2e_envvars, cint_g2e_index_xyz
   public :: cint_g0_2e_2d, cint_g0_2e
   public :: cint_g0_lj2d_4d, cint_g0_kj2d_4d, cint_g0_il2d_4d, cint_g0_ik2d_4d
   public :: cint_all_2e_optimizer
   public :: cint_all_3c2e_optimizer
   public :: cint_all_2c2e_optimizer

   real(dp), parameter :: SQRTPI = 1.7724538509055160272981674833411451_dp
   real(dp), parameter :: PI     = 3.1415926535897932384626433832795029_dp
   real(dp), parameter :: EXPCUTOFF_SR = 40.0_dp

contains

   subroutine cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      type(cint_env_vars), intent(inout) :: envs
      integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
      ! no INTENT(IN) on the three tables: envs points at them, and F2018
      ! 8.5.10 forbids an INTENT(IN) dummy as a pointer target.  They are
      ! never written -- see the note on cint_env_vars.
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)

      integer  :: i_sh, j_sh, k_sh, l_sh, dli, dlj, dlk, dll, ip
      integer  :: rys_order, nrys_roots
      logical  :: ibase, kbase
      real(dp) :: omega

      envs%natm = natm
      envs%nbas = nbas
      envs%atm(0:) => atm
      envs%bas(0:) => bas
      envs%env(0:) => env
      envs%shls(0:size(shls)-1) = shls
      envs%ng(0:min(7,size(ng)-1)) = ng(0:min(7,size(ng)-1))
      envs%omega = env(PTR_RANGE_OMEGA)

      i_sh = shls(0); j_sh = shls(1); k_sh = shls(2); l_sh = shls(3)
      envs%i_l = bas_of(envs, ANG_OF, i_sh)
      envs%j_l = bas_of(envs, ANG_OF, j_sh)
      envs%k_l = bas_of(envs, ANG_OF, k_sh)
      envs%l_l = bas_of(envs, ANG_OF, l_sh)
      envs%x_ctr(0) = bas_of(envs, NCTR_OF, i_sh)
      envs%x_ctr(1) = bas_of(envs, NCTR_OF, j_sh)
      envs%x_ctr(2) = bas_of(envs, NCTR_OF, k_sh)
      envs%x_ctr(3) = bas_of(envs, NCTR_OF, l_sh)
      ! L SHELLS.  ANG_OF already gave i_l = 1 for one, so every ceiling,
      ! stride and Rys order below is the p sub-block's and needs no change
      ! -- the g tensor at li_ceil = 1 holds both components already, index 0
      ! being the s function and the x/y/z blocks at index 1 the p.  What
      ! does change is the component count: four, not three.
      envs%sp_mask = 0
      if (cint_bas_is_sp(i_sh, bas)) envs%sp_mask = ibset(envs%sp_mask, 0)
      if (cint_bas_is_sp(j_sh, bas)) envs%sp_mask = ibset(envs%sp_mask, 1)
      if (cint_bas_is_sp(k_sh, bas)) envs%sp_mask = ibset(envs%sp_mask, 2)
      if (cint_bas_is_sp(l_sh, bas)) envs%sp_mask = ibset(envs%sp_mask, 3)
      envs%nfi = (envs%i_l+1)*(envs%i_l+2)/2
      envs%nfj = (envs%j_l+1)*(envs%j_l+2)/2
      envs%nfk = (envs%k_l+1)*(envs%k_l+2)/2
      envs%nfl = (envs%l_l+1)*(envs%l_l+2)/2
      if (envs%sp_mask /= 0) then
         if (btest(envs%sp_mask, 0)) envs%nfi = NF_SP
         if (btest(envs%sp_mask, 1)) envs%nfj = NF_SP
         if (btest(envs%sp_mask, 2)) envs%nfk = NF_SP
         if (btest(envs%sp_mask, 3)) envs%nfl = NF_SP
      end if
      ! note the order: i, k, l, j -- the g array is laid out (i,k,l,j)
      envs%nf = envs%nfi * envs%nfk * envs%nfl * envs%nfj

      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, i_sh))
      envs%ri(0) = env(ip); envs%ri(1) = env(ip+1); envs%ri(2) = env(ip+2)
      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, j_sh))
      envs%rj(0) = env(ip); envs%rj(1) = env(ip+1); envs%rj(2) = env(ip+2)
      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, k_sh))
      envs%rk(0) = env(ip); envs%rk(1) = env(ip+1); envs%rk(2) = env(ip+2)
      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, l_sh))
      envs%rl(0) = env(ip); envs%rl(1) = env(ip+1); envs%rl(2) = env(ip+2)

      envs%common_factor = (PI*PI*PI) * 2.0_dp / SQRTPI &
         * sp_split_fac(envs, 0, envs%i_l) * sp_split_fac(envs, 1, envs%j_l) &
         * sp_split_fac(envs, 2, envs%k_l) * sp_split_fac(envs, 3, envs%l_l)
      if (env(PTR_EXPCUTOFF) == 0.0_dp) then
         envs%expcutoff = EXPCUTOFF
      else
         ! +1 to ensure accuracy, as the C notes
         envs%expcutoff = max(MIN_EXPCUTOFF, env(PTR_EXPCUTOFF)) + 1.0_dp
      end if

      envs%gbits = ng(GSHIFT)
      envs%ncomp_e1 = ng(POS_E1)
      envs%ncomp_e2 = ng(POS_E2)
      envs%ncomp_tensor = ng(TENSOR)

      envs%li_ceil = envs%i_l + ng(IINC)
      envs%lj_ceil = envs%j_l + ng(JINC)
      envs%lk_ceil = envs%k_l + ng(KINC)
      envs%ll_ceil = envs%l_l + ng(LINC)
      rys_order = (envs%li_ceil + envs%lj_ceil + envs%lk_ceil + envs%ll_ceil)/2 + 1
      nrys_roots = rys_order
      omega = env(PTR_RANGE_OMEGA)
      ! the short-range kernel needs twice the roots at low order
      if (omega < 0.0_dp .and. rys_order <= 3) nrys_roots = nrys_roots * 2
      envs%rys_order = rys_order
      envs%nrys_roots = nrys_roots

      ibase = envs%li_ceil > envs%lj_ceil
      kbase = envs%lk_ceil > envs%ll_ceil
      if (kbase) then
         dlk = envs%lk_ceil + envs%ll_ceil + 1
         dll = envs%ll_ceil + 1
      else
         dlk = envs%lk_ceil + 1
         dll = envs%lk_ceil + envs%ll_ceil + 1
      end if
      if (ibase) then
         dli = envs%li_ceil + envs%lj_ceil + 1
         dlj = envs%lj_ceil + 1
      else
         dli = envs%li_ceil + 1
         dlj = envs%li_ceil + envs%lj_ceil + 1
      end if
      envs%g_stride_i = nrys_roots
      envs%g_stride_k = nrys_roots * dli
      envs%g_stride_l = nrys_roots * dli * dlk
      envs%g_stride_j = nrys_roots * dli * dlk * dll
      envs%g_size     = nrys_roots * dli * dlk * dll * dlj

      if (kbase) then
         envs%g2d_klmax = envs%g_stride_k
         envs%rx_in_rklrx = envs%rk
         envs%rkrl = envs%rk - envs%rl
      else
         envs%g2d_klmax = envs%g_stride_l
         envs%rx_in_rklrx = envs%rl
         envs%rkrl = envs%rl - envs%rk
      end if

      if (ibase) then
         envs%g2d_ijmax = envs%g_stride_i
         envs%rx_in_rijrx = envs%ri
         envs%rirj = envs%ri - envs%rj
      else
         envs%g2d_ijmax = envs%g_stride_j
         envs%rx_in_rijrx = envs%rj
         envs%rirj = envs%rj - envs%ri
      end if

      ! Up to two Rys roots the transfer is done by a table of unrolled
      ! kernels; above that by the general 2D-to-4D routines, chosen by which
      ! centre each pair is built along.
      if (rys_order <= 2) then
         envs%f_g0_2d4d = G2D4D_UNROLLED
         if (rys_order /= nrys_roots) envs%f_g0_2d4d = G2D4D_SR_UNROLLED
      else if (kbase) then
         if (ibase) then
            envs%f_g0_2d4d = G2D4D_IK
         else
            envs%f_g0_2d4d = G2D4D_KJ
         end if
      else
         if (ibase) then
            envs%f_g0_2d4d = G2D4D_IL
         else
            envs%f_g0_2d4d = G2D4D_LJ
         end if
      end if
   end subroutine cint_init_int2e_envvars

   ! Three-centre two-electron: (ij|k).  The C's trick, kept here, is to load
   ! the auxiliary shell into the *l* slot and leave k at zero, so the whole
   ! 2D-to-4D machinery of the four-centre case is reused unchanged -- hence
   ! `ll_ceil = k_l + ng(KINC)` below, which reads like a typo and is not.
   subroutine cint_init_int3c2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      type(cint_env_vars), intent(inout) :: envs
      integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
      ! no INTENT(IN) on the three tables: envs points at them, and F2018
      ! 8.5.10 forbids an INTENT(IN) dummy as a pointer target.  They are
      ! never written -- see the note on cint_env_vars.
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)

      integer  :: i_sh, j_sh, k_sh, dli, dlj, dlk, ip
      integer  :: rys_order, nrys_roots
      logical  :: ibase
      real(dp) :: omega

      envs%natm = natm
      envs%nbas = nbas
      envs%atm(0:) => atm
      envs%bas(0:) => bas
      envs%env(0:) => env
      envs%shls(0:size(shls)-1) = shls
      envs%ng(0:min(7,size(ng)-1)) = ng(0:min(7,size(ng)-1))
      envs%omega = env(PTR_RANGE_OMEGA)

      i_sh = shls(0); j_sh = shls(1); k_sh = shls(2)
      envs%i_l = bas_of(envs, ANG_OF, i_sh)
      envs%j_l = bas_of(envs, ANG_OF, j_sh)
      envs%k_l = bas_of(envs, ANG_OF, k_sh)
      envs%l_l = 0
      envs%x_ctr(0) = bas_of(envs, NCTR_OF, i_sh)
      envs%x_ctr(1) = bas_of(envs, NCTR_OF, j_sh)
      envs%x_ctr(2) = bas_of(envs, NCTR_OF, k_sh)
      envs%x_ctr(3) = 1
      envs%nfi = (envs%i_l+1)*(envs%i_l+2)/2
      envs%nfj = (envs%j_l+1)*(envs%j_l+2)/2
      envs%nfk = (envs%k_l+1)*(envs%k_l+2)/2
      envs%nfl = 1
      envs%nf = envs%nfi * envs%nfk * envs%nfj

      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, i_sh))
      envs%ri(0) = env(ip); envs%ri(1) = env(ip+1); envs%ri(2) = env(ip+2)
      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, j_sh))
      envs%rj(0) = env(ip); envs%rj(1) = env(ip+1); envs%rj(2) = env(ip+2)
      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, k_sh))
      envs%rk(0) = env(ip); envs%rk(1) = env(ip+1); envs%rk(2) = env(ip+2)

      envs%common_factor = (PI*PI*PI) * 2.0_dp / SQRTPI &
         * cint_common_fac_sp(envs%i_l) * cint_common_fac_sp(envs%j_l) &
         * cint_common_fac_sp(envs%k_l)
      ! note: no +1 here, unlike the four-centre case
      if (env(PTR_EXPCUTOFF) == 0.0_dp) then
         envs%expcutoff = EXPCUTOFF
      else
         envs%expcutoff = max(MIN_EXPCUTOFF, env(PTR_EXPCUTOFF))
      end if

      envs%gbits = ng(GSHIFT)
      envs%ncomp_e1 = ng(POS_E1)
      envs%ncomp_e2 = ng(POS_E2)
      envs%ncomp_tensor = ng(TENSOR)

      envs%li_ceil = envs%i_l + ng(IINC)
      envs%lj_ceil = envs%j_l + ng(JINC)
      envs%lk_ceil = 0
      envs%ll_ceil = envs%k_l + ng(KINC)
      rys_order = (envs%li_ceil + envs%lj_ceil + envs%ll_ceil)/2 + 1
      nrys_roots = rys_order
      omega = env(PTR_RANGE_OMEGA)
      if (omega < 0.0_dp .and. rys_order <= 3) nrys_roots = nrys_roots * 2
      envs%rys_order = rys_order
      envs%nrys_roots = nrys_roots

      ibase = envs%li_ceil > envs%lj_ceil
      if (ibase) then
         dli = envs%li_ceil + envs%lj_ceil + 1
         dlj = envs%lj_ceil + 1
      else
         dli = envs%li_ceil + 1
         dlj = envs%li_ceil + envs%lj_ceil + 1
      end if
      dlk = envs%ll_ceil + 1

      envs%g_stride_i = nrys_roots
      envs%g_stride_k = nrys_roots * dli
      envs%g_stride_l = nrys_roots * dli
      envs%g_stride_j = nrys_roots * dli * dlk
      envs%g_size     = nrys_roots * dli * dlk * dlj

      envs%al = 0.0_dp
      envs%rkl = envs%rk
      envs%g2d_klmax = envs%g_stride_k
      ! in g0_2d rklrx = rkl - rx = 0, so rkl is the reference point
      envs%rkrl = envs%rk
      envs%rx_in_rklrx = envs%rk

      if (ibase) then
         envs%g2d_ijmax = envs%g_stride_i
         envs%rx_in_rijrx = envs%ri
         envs%rirj = envs%ri - envs%rj
      else
         envs%g2d_ijmax = envs%g_stride_j
         envs%rx_in_rijrx = envs%rj
         envs%rirj = envs%rj - envs%ri
      end if

      if (rys_order <= 2) then
         envs%f_g0_2d4d = G2D4D_UNROLLED
         if (rys_order /= nrys_roots) envs%f_g0_2d4d = G2D4D_SR_UNROLLED
      else if (ibase) then
         envs%f_g0_2d4d = G2D4D_IL
      else
         envs%f_g0_2d4d = G2D4D_LJ
      end if
   end subroutine cint_init_int3c2e_envvars

   ! Two-centre two-electron: (i|k).  Only two indices survive, so the 2D
   ! recursion output needs no transfer at all.
   subroutine cint_init_int2c2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      type(cint_env_vars), intent(inout) :: envs
      integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
      ! no INTENT(IN) on the three tables: envs points at them, and F2018
      ! 8.5.10 forbids an INTENT(IN) dummy as a pointer target.  They are
      ! never written -- see the note on cint_env_vars.
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)

      integer  :: i_sh, k_sh, dli, dlk, ip
      integer  :: rys_order, nrys_roots
      real(dp) :: omega

      envs%natm = natm
      envs%nbas = nbas
      envs%atm(0:) => atm
      envs%bas(0:) => bas
      envs%env(0:) => env
      envs%shls(0:size(shls)-1) = shls
      envs%ng(0:min(7,size(ng)-1)) = ng(0:min(7,size(ng)-1))
      envs%omega = env(PTR_RANGE_OMEGA)

      i_sh = shls(0); k_sh = shls(1)
      envs%i_l = bas_of(envs, ANG_OF, i_sh)
      envs%j_l = 0
      envs%k_l = bas_of(envs, ANG_OF, k_sh)
      envs%l_l = 0
      envs%x_ctr(0) = bas_of(envs, NCTR_OF, i_sh)
      envs%x_ctr(1) = bas_of(envs, NCTR_OF, k_sh)
      envs%x_ctr(2) = 1
      envs%x_ctr(3) = 1
      envs%nfi = (envs%i_l+1)*(envs%i_l+2)/2
      envs%nfj = 1
      envs%nfk = (envs%k_l+1)*(envs%k_l+2)/2
      envs%nfl = 1
      envs%nf = envs%nfi * envs%nfk

      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, i_sh))
      envs%ri(0) = env(ip); envs%ri(1) = env(ip+1); envs%ri(2) = env(ip+2)
      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, k_sh))
      envs%rk(0) = env(ip); envs%rk(1) = env(ip+1); envs%rk(2) = env(ip+2)

      envs%common_factor = (PI*PI*PI) * 2.0_dp / SQRTPI &
         * cint_common_fac_sp(envs%i_l) * cint_common_fac_sp(envs%k_l)
      if (env(PTR_EXPCUTOFF) == 0.0_dp) then
         envs%expcutoff = EXPCUTOFF
      else
         envs%expcutoff = max(MIN_EXPCUTOFF, env(PTR_EXPCUTOFF))
      end if

      envs%gbits = ng(GSHIFT)
      envs%ncomp_e1 = ng(POS_E1)
      envs%ncomp_e2 = ng(POS_E2)
      envs%ncomp_tensor = ng(TENSOR)

      envs%li_ceil = envs%i_l + ng(IINC)
      envs%lj_ceil = 0
      envs%lk_ceil = envs%k_l + ng(KINC)
      envs%ll_ceil = 0
      rys_order = (envs%li_ceil + envs%lk_ceil)/2 + 1
      nrys_roots = rys_order
      omega = env(PTR_RANGE_OMEGA)
      if (omega < 0.0_dp .and. rys_order <= 3) nrys_roots = nrys_roots * 2
      envs%rys_order = rys_order
      envs%nrys_roots = nrys_roots

      dli = envs%li_ceil + 1
      dlk = envs%lk_ceil + 1
      envs%g_stride_i = nrys_roots
      envs%g_stride_k = nrys_roots * dli
      envs%g_stride_l = envs%g_stride_k
      envs%g_size     = nrys_roots * dli * dlk

      envs%aj = 0.0_dp
      envs%al = 0.0_dp
      envs%rij = envs%ri
      envs%rkl = envs%rk
      envs%g2d_ijmax = envs%g_stride_i
      envs%g2d_klmax = envs%g_stride_k
      envs%rkrl = envs%rk
      envs%rirj = envs%ri
      envs%rx_in_rklrx = envs%rk
      envs%rx_in_rijrx = envs%ri

      if (rys_order <= 2) then
         envs%f_g0_2d4d = G2D4D_UNROLLED
         if (rys_order /= nrys_roots) envs%f_g0_2d4d = G2D4D_SR_UNROLLED
      else
         envs%f_g0_2d4d = G2D4D_NONE
      end if

      ! j_l, nfj and g_stride_j are set from the k values because c2s_sph_1e
      ! and CINTg1e_index_xyz -- which int2c2e reuses -- address the second
      ! index as if it were j.
      envs%j_l = envs%k_l
      envs%nfj = envs%nfk
      envs%g_stride_j = envs%g_stride_k
   end subroutine cint_init_int2c2e_envvars

   ! The shell-normalisation factor common_factor carries for one index.
   !
   ! An L shell's two sub-blocks want different ones -- 1/sqrt(4pi) for the
   ! s, sqrt(3/4pi) for the p -- and common_factor is a single scalar for the
   ! whole quartet, so it cannot carry either.  It carries 1 instead and the
   ! per-component factor rides on the contraction coefficient, in
   ! cint_prim_to_ctr_sp_0/_1, which is the only other place a per-component
   ! scalar can be applied without touching the g tensor.
   pure function sp_split_fac(envs, x, l) result(f)
      type(cint_env_vars), intent(in) :: envs
      integer,             intent(in) :: x, l
      real(dp) :: f
      if (btest(envs%sp_mask, x)) then
         f = 1.0_dp
      else
         f = cint_common_fac_sp(l)
      end if
   end function sp_split_fac

   ! Cartesian component powers for one index of a quartet, L shell or not.
   pure subroutine sp_cart_comp(nx, ny, nz, envs, x, l)
      integer, intent(out) :: nx(0:), ny(0:), nz(0:)
      type(cint_env_vars), intent(in) :: envs
      integer, intent(in) :: x, l
      if (btest(envs%sp_mask, x)) then
         call cint_cart_comp_sp(nx, ny, nz)
      else
         call cint_cart_comp(nx, ny, nz, l)
      end if
   end subroutine sp_cart_comp

   ! The index map when at least one of the four shells is an L shell.
   !
   ! Structurally the routine below with the unrolling taken out -- and a
   ! separate routine for exactly that reason.  The i_l select there writes
   ! the s, p and d offsets as immediates because that is what every basis
   ! set is made of; threading a fifth case through it would put a test in
   ! front of the case that matters.  An L shell's four components are the
   ! l = 1 map with (0,0,0) prepended, which the general form already
   ! produces once cint_cart_comp_sp supplies the table.
   subroutine g2e_index_xyz_sp(idx, envs)
      integer, intent(out) :: idx(0:)
      type(cint_env_vars), intent(in) :: envs

      integer :: i, j, k, l, n, di, dj, dk, dl, ofx, ofy, ofz
      integer :: ofkx, ofky, ofkz, oflx, ofly, oflz, ofjx, ofjy, ofjz
      integer :: i_nx(0:CART_MAX-1), i_ny(0:CART_MAX-1), i_nz(0:CART_MAX-1)
      integer :: j_nx(0:CART_MAX-1), j_ny(0:CART_MAX-1), j_nz(0:CART_MAX-1)
      integer :: k_nx(0:CART_MAX-1), k_ny(0:CART_MAX-1), k_nz(0:CART_MAX-1)
      integer :: l_nx(0:CART_MAX-1), l_ny(0:CART_MAX-1), l_nz(0:CART_MAX-1)

      di = envs%g_stride_i
      dk = envs%g_stride_k
      dl = envs%g_stride_l
      dj = envs%g_stride_j
      call sp_cart_comp(i_nx, i_ny, i_nz, envs, 0, envs%i_l)
      call sp_cart_comp(j_nx, j_ny, j_nz, envs, 1, envs%j_l)
      call sp_cart_comp(k_nx, k_ny, k_nz, envs, 2, envs%k_l)
      call sp_cart_comp(l_nx, l_ny, l_nz, envs, 3, envs%l_l)

      ofx = 0
      ofy = envs%g_size
      ofz = envs%g_size * 2
      n = 0
      do j = 0, envs%nfj - 1
         ofjx = ofx + dj * j_nx(j)
         ofjy = ofy + dj * j_ny(j)
         ofjz = ofz + dj * j_nz(j)
         do l = 0, envs%nfl - 1
            oflx = ofjx + dl * l_nx(l)
            ofly = ofjy + dl * l_ny(l)
            oflz = ofjz + dl * l_nz(l)
            do k = 0, envs%nfk - 1
               ofkx = oflx + dk * k_nx(k)
               ofky = ofly + dk * k_ny(k)
               ofkz = oflz + dk * k_nz(k)
               do i = 0, envs%nfi - 1
                  idx(n+0) = ofkx + di * i_nx(i)
                  idx(n+1) = ofky + di * i_ny(i)
                  idx(n+2) = ofkz + di * i_nz(i)
                  n = n + 3
               end do
            end do
         end do
      end do
   end subroutine g2e_index_xyz_sp

   ! For each (i,k,l,j) Cartesian component, the three offsets into g.
   subroutine cint_g2e_index_xyz(idx, envs)
      integer, intent(out) :: idx(0:)
      type(cint_env_vars), intent(in) :: envs

      integer :: i, j, k, l, n, di, dj, dk, dl, ofx, ofy, ofz
      integer :: ofkx, ofky, ofkz, oflx, ofly, oflz, ofjx, ofjy, ofjz
      integer :: i_nx(0:CART_MAX-1), i_ny(0:CART_MAX-1), i_nz(0:CART_MAX-1)
      integer :: j_nx(0:CART_MAX-1), j_ny(0:CART_MAX-1), j_nz(0:CART_MAX-1)
      integer :: k_nx(0:CART_MAX-1), k_ny(0:CART_MAX-1), k_nz(0:CART_MAX-1)
      integer :: l_nx(0:CART_MAX-1), l_ny(0:CART_MAX-1), l_nz(0:CART_MAX-1)

      if (envs%sp_mask /= 0) then
         call g2e_index_xyz_sp(idx, envs)
         return
      end if

      di = envs%g_stride_i
      dk = envs%g_stride_k
      dl = envs%g_stride_l
      dj = envs%g_stride_j
      call cint_cart_comp(i_nx, i_ny, i_nz, envs%i_l)
      call cint_cart_comp(j_nx, j_ny, j_nz, envs%j_l)
      call cint_cart_comp(k_nx, k_ny, k_nz, envs%k_l)
      call cint_cart_comp(l_nx, l_ny, l_nz, envs%l_l)

      ofx = 0
      ofy = envs%g_size
      ofz = envs%g_size * 2
      n = 0
      do j = 0, envs%nfj - 1
         ofjx = ofx + dj * j_nx(j)
         ofjy = ofy + dj * j_ny(j)
         ofjz = ofz + dj * j_nz(j)
         do l = 0, envs%nfl - 1
            oflx = ofjx + dl * l_nx(l)
            ofly = ofjy + dl * l_ny(l)
            oflz = ofjz + dl * l_nz(l)
            do k = 0, envs%nfk - 1
               ofkx = oflx + dk * k_nx(k)
               ofky = ofly + dk * k_ny(k)
               ofkz = oflz + dk * k_nz(k)
               ! UNROLLED ON i_l, as the C is (src/g2e.c).  s, p and d are
               ! most of every basis set, and for those the C writes the
               ! offsets as immediates instead of reading i_nx, i_ny, i_nz
               ! back out of the tables it just built.  Isolated, the two
               ! were 1.54x apart at p and 1.75x at d, and level at f and g
               ! -- which is where the C falls through to the same loop.
               select case (envs%i_l)
               case (0)
                  idx(n+0) = ofkx
                  idx(n+1) = ofky
                  idx(n+2) = ofkz
                  n = n + 3
               case (1)
                  idx(n+0) = ofkx + di
                  idx(n+1) = ofky
                  idx(n+2) = ofkz
                  idx(n+3) = ofkx
                  idx(n+4) = ofky + di
                  idx(n+5) = ofkz
                  idx(n+6) = ofkx
                  idx(n+7) = ofky
                  idx(n+8) = ofkz + di
                  n = n + 9
               case (2)
                  idx(n+0)  = ofkx + di*2
                  idx(n+1)  = ofky
                  idx(n+2)  = ofkz
                  idx(n+3)  = ofkx + di
                  idx(n+4)  = ofky + di
                  idx(n+5)  = ofkz
                  idx(n+6)  = ofkx + di
                  idx(n+7)  = ofky
                  idx(n+8)  = ofkz + di
                  idx(n+9)  = ofkx
                  idx(n+10) = ofky + di*2
                  idx(n+11) = ofkz
                  idx(n+12) = ofkx
                  idx(n+13) = ofky + di
                  idx(n+14) = ofkz + di
                  idx(n+15) = ofkx
                  idx(n+16) = ofky
                  idx(n+17) = ofkz + di*2
                  n = n + 18
               case default
                  do i = 0, envs%nfi - 1
                     idx(n+0) = ofkx + di * i_nx(i)
                     idx(n+1) = ofky + di * i_ny(i)
                     idx(n+2) = ofkz + di * i_nz(i)
                     n = n + 3
                  end do
               end select
            end do
         end do
      end do
   end subroutine cint_g2e_index_xyz

   ! The Rys 2D recursion.  Builds g(m,n) for m along the kl pair and n along
   ! the ij pair, one root at a time, using the b/c coefficients.
   subroutine cint_g0_2e_2d(g, bc, envs)
      real(dp), intent(inout), target :: g(0:)
      type(rys_2e_t), intent(in) :: bc
      type(cint_env_vars), intent(in) :: envs

      integer  :: nroots, nmax, mmax, dm, dn, i, j, m, n, off, gx, gy, gz
      real(dp) :: s0x, s1x, s2x, s0y, s1y, s2y, s0z, s1z, s2z
      real(dp) :: c00x, c00y, c00z, c0px, c0py, c0pz, b10, b01, b00
      ! THREE ROWS, BOUND ONCE.  The C's DEF_GXYZ hands this routine three
      ! separate pointers into g, so gx[j] is one indexed load; g(gx+j) is
      ! that plus an add, on every one of the recursion's accesses.  It was
      ! the routine's largest single line at 17.2 million instructions.
      real(dp), pointer, contiguous :: gxp(:), gyp(:), gzp(:)

      nroots = envs%nrys_roots
      nmax = envs%li_ceil + envs%lj_ceil
      mmax = envs%lk_ceil + envs%ll_ceil
      dm = envs%g2d_klmax
      dn = envs%g2d_ijmax
      gx = 0; gy = envs%g_size; gz = envs%g_size * 2
      gxp(0:) => g(gx:); gyp(0:) => g(gy:); gzp(0:) => g(gz:)

      do i = 0, nroots - 1
         gxp(i) = 1.0_dp
         gyp(i) = 1.0_dp
         ! gzp(i) is the Rys weight, already there
      end do

      do i = 0, nroots - 1
         c00x = bc%c00x(i); c00y = bc%c00y(i); c00z = bc%c00z(i)
         c0px = bc%c0px(i); c0py = bc%c0py(i); c0pz = bc%c0pz(i)
         b10 = bc%b10(i); b01 = bc%b01(i); b00 = bc%b00(i)

         if (nmax > 0) then
            ! gx(0,n+1) = c00*gx(0,n) + n*b10*gx(0,n-1)
            s0x = gxp(i); s0y = gyp(i); s0z = gzp(i)
            s1x = c00x * s0x; s1y = c00y * s0y; s1z = c00z * s0z
            gxp(i+dn) = s1x; gyp(i+dn) = s1y; gzp(i+dn) = s1z
            do n = 1, nmax - 1
               s2x = c00x * s1x + n * b10 * s0x
               s2y = c00y * s1y + n * b10 * s0y
               s2z = c00z * s1z + n * b10 * s0z
               gxp(i+(n+1)*dn) = s2x
               gyp(i+(n+1)*dn) = s2y
               gzp(i+(n+1)*dn) = s2z
               s0x = s1x; s0y = s1y; s0z = s1z
               s1x = s2x; s1y = s2y; s1z = s2z
            end do
         end if

         if (mmax > 0) then
            ! gx(m+1,0) = c0p*gx(m,0) + m*b01*gx(m-1,0)
            s0x = gxp(i); s0y = gyp(i); s0z = gzp(i)
            s1x = c0px * s0x; s1y = c0py * s0y; s1z = c0pz * s0z
            gxp(i+dm) = s1x; gyp(i+dm) = s1y; gzp(i+dm) = s1z
            do m = 1, mmax - 1
               s2x = c0px * s1x + m * b01 * s0x
               s2y = c0py * s1y + m * b01 * s0y
               s2z = c0pz * s1z + m * b01 * s0z
               gxp(i+(m+1)*dm) = s2x
               gyp(i+(m+1)*dm) = s2y
               gzp(i+(m+1)*dm) = s2z
               s0x = s1x; s0y = s1y; s0z = s1z
               s1x = s2x; s1y = s2y; s1z = s2z
            end do

            if (nmax > 0) then
               ! gx(m+1,1) = c0p*gx(m,1) + m*b01*gx(m-1,1) + b00*gx(m,0)
               s0x = gxp(i+dn); s0y = gyp(i+dn); s0z = gzp(i+dn)
               s1x = c0px * s0x + b00 * gxp(i)
               s1y = c0py * s0y + b00 * gyp(i)
               s1z = c0pz * s0z + b00 * gzp(i)
               gxp(i+dn+dm) = s1x; gyp(i+dn+dm) = s1y; gzp(i+dn+dm) = s1z
               do m = 1, mmax - 1
                  s2x = c0px*s1x + m*b01*s0x + b00*gxp(i+m*dm)
                  s2y = c0py*s1y + m*b01*s0y + b00*gyp(i+m*dm)
                  s2z = c0pz*s1z + m*b01*s0z + b00*gzp(i+m*dm)
                  gxp(i+dn+(m+1)*dm) = s2x
                  gyp(i+dn+(m+1)*dm) = s2y
                  gzp(i+dn+(m+1)*dm) = s2z
                  s0x = s1x; s0y = s1y; s0z = s1z
                  s1x = s2x; s1y = s2y; s1z = s2z
               end do
            end if
         end if

         ! gx(m,n+1) = c00*gx(m,n) + n*b10*gx(m,n-1) + m*b00*gx(m-1,n)
         do m = 1, mmax
            off = m * dm
            j = off + i
            s0x = gxp(j); s0y = gyp(j); s0z = gzp(j)
            s1x = gxp(j+dn); s1y = gyp(j+dn); s1z = gzp(j+dn)
            do n = 1, nmax - 1
               s2x = c00x*s1x + n*b10*s0x + m*b00*gxp(j+n*dn-dm)
               s2y = c00y*s1y + n*b10*s0y + m*b00*gyp(j+n*dn-dm)
               s2z = c00z*s1z + n*b10*s0z + m*b00*gzp(j+n*dn-dm)
               gxp(j+(n+1)*dn) = s2x
               gyp(j+(n+1)*dn) = s2y
               gzp(j+(n+1)*dn) = s2z
               s0x = s1x; s0y = s1y; s0z = s1z
               s1x = s2x; s1y = s2y; s1z = s2z
            end do
         end do
      end do
   end subroutine cint_g0_2e_2d

   ! The four transfers.  Each shifts angular momentum between the two centres
   ! of a pair, using g(i,..) = rirj*g(i-1,..) + g(i-1,..,j+1) and its kl
   ! analogue.  The C writes them with p1/p2 pointers offset by -di etc.;
   ! those become explicit index shifts here.

   subroutine cint_g0_lj2d_4d(g, envs)
      real(dp), intent(inout), target :: g(0:)
      ! Three rows bound once, as the C's DEF_GXYZ does; see cint_g0_2e_2d.
      real(dp), pointer, contiguous :: gxp(:), gyp(:), gzp(:)
      type(cint_env_vars), intent(in) :: envs
      integer  :: li, lk, lj, nmax, mmax, nroots, i, j, k, l, ptr, n
      integer  :: di, dk, dl, dj, gx, gy, gz
      real(dp) :: rx, ry, rz

      li = envs%li_ceil; lk = envs%lk_ceil
      if (li == 0 .and. lk == 0) return
      nmax = envs%li_ceil + envs%lj_ceil
      mmax = envs%lk_ceil + envs%ll_ceil
      lj = envs%lj_ceil
      nroots = envs%nrys_roots
      di = envs%g_stride_i; dk = envs%g_stride_k
      dl = envs%g_stride_l; dj = envs%g_stride_j
      gx = 0; gy = envs%g_size; gz = envs%g_size*2
      gxp(0:) => g(gx:); gyp(0:) => g(gy:); gzp(0:) => g(gz:)

      rx = envs%rirj(0); ry = envs%rirj(1); rz = envs%rirj(2)
      do i = 1, li
         do j = 0, nmax - i
            do l = 0, mmax
               ptr = j*dj + l*dl + i*di
               do n = ptr, ptr + nroots - 1
                  gxp(n) = rx * gxp(n-di) + gxp(n-di+dj)
                  gyp(n) = ry * gyp(n-di) + gyp(n-di+dj)
                  gzp(n) = rz * gzp(n-di) + gzp(n-di+dj)
               end do
            end do
         end do
      end do

      rx = envs%rkrl(0); ry = envs%rkrl(1); rz = envs%rkrl(2)
      do j = 0, lj
         do k = 1, lk
            do l = 0, mmax - k
               ptr = j*dj + l*dl + k*dk
               do n = ptr, ptr + dk - 1
                  gxp(n) = rx * gxp(n-dk) + gxp(n-dk+dl)
                  gyp(n) = ry * gyp(n-dk) + gyp(n-dk+dl)
                  gzp(n) = rz * gzp(n-dk) + gzp(n-dk+dl)
               end do
            end do
         end do
      end do
   end subroutine cint_g0_lj2d_4d

   subroutine cint_g0_kj2d_4d(g, envs)
      real(dp), intent(inout), target :: g(0:)
      ! Three rows bound once, as the C's DEF_GXYZ does; see cint_g0_2e_2d.
      real(dp), pointer, contiguous :: gxp(:), gyp(:), gzp(:)
      type(cint_env_vars), intent(in) :: envs
      integer  :: li, ll, lj, nmax, mmax, nroots, i, j, k, l, ptr, n
      integer  :: di, dk, dl, dj, gx, gy, gz
      real(dp) :: rx, ry, rz

      li = envs%li_ceil; ll = envs%ll_ceil
      if (li == 0 .and. ll == 0) return
      nmax = envs%li_ceil + envs%lj_ceil
      mmax = envs%lk_ceil + envs%ll_ceil
      lj = envs%lj_ceil
      nroots = envs%nrys_roots
      di = envs%g_stride_i; dk = envs%g_stride_k
      dl = envs%g_stride_l; dj = envs%g_stride_j
      gx = 0; gy = envs%g_size; gz = envs%g_size*2
      gxp(0:) => g(gx:); gyp(0:) => g(gy:); gzp(0:) => g(gz:)

      rx = envs%rirj(0); ry = envs%rirj(1); rz = envs%rirj(2)
      do i = 1, li
         do j = 0, nmax - i
            do k = 0, mmax
               ptr = j*dj + k*dk + i*di
               do n = ptr, ptr + nroots - 1
                  gxp(n) = rx * gxp(n-di) + gxp(n-di+dj)
                  gyp(n) = ry * gyp(n-di) + gyp(n-di+dj)
                  gzp(n) = rz * gzp(n-di) + gzp(n-di+dj)
               end do
            end do
         end do
      end do

      rx = envs%rkrl(0); ry = envs%rkrl(1); rz = envs%rkrl(2)
      do j = 0, lj
         do l = 1, ll
            do k = 0, mmax - l
               ptr = j*dj + l*dl + k*dk
               do n = ptr, ptr + dk - 1
                  gxp(n) = rx * gxp(n-dl) + gxp(n-dl+dk)
                  gyp(n) = ry * gyp(n-dl) + gyp(n-dl+dk)
                  gzp(n) = rz * gzp(n-dl) + gzp(n-dl+dk)
               end do
            end do
         end do
      end do
   end subroutine cint_g0_kj2d_4d

   subroutine cint_g0_il2d_4d(g, envs)
      real(dp), intent(inout), target :: g(0:)
      ! Three rows bound once, as the C's DEF_GXYZ does; see cint_g0_2e_2d.
      real(dp), pointer, contiguous :: gxp(:), gyp(:), gzp(:)
      type(cint_env_vars), intent(in) :: envs
      integer  :: lk, lj, ll, nmax, mmax, nroots, i, j, k, l, ptr, n
      integer  :: di, dk, dl, dj, gx, gy, gz
      real(dp) :: rx, ry, rz

      lk = envs%lk_ceil; lj = envs%lj_ceil
      if (lj == 0 .and. lk == 0) return
      nmax = envs%li_ceil + envs%lj_ceil
      mmax = envs%lk_ceil + envs%ll_ceil
      ll = envs%ll_ceil
      nroots = envs%nrys_roots
      di = envs%g_stride_i; dk = envs%g_stride_k
      dl = envs%g_stride_l; dj = envs%g_stride_j
      gx = 0; gy = envs%g_size; gz = envs%g_size*2
      gxp(0:) => g(gx:); gyp(0:) => g(gy:); gzp(0:) => g(gz:)

      rx = envs%rkrl(0); ry = envs%rkrl(1); rz = envs%rkrl(2)
      do k = 1, lk
         do l = 0, mmax - k
            do i = 0, nmax
               ptr = l*dl + k*dk + i*di
               do n = ptr, ptr + nroots - 1
                  gxp(n) = rx * gxp(n-dk) + gxp(n-dk+dl)
                  gyp(n) = ry * gyp(n-dk) + gyp(n-dk+dl)
                  gzp(n) = rz * gzp(n-dk) + gzp(n-dk+dl)
               end do
            end do
         end do
      end do

      rx = envs%rirj(0); ry = envs%rirj(1); rz = envs%rirj(2)
      do j = 1, lj
         do l = 0, ll
            do k = 0, lk
               ptr = j*dj + l*dl + k*dk
               do n = ptr, ptr + dk - di*j - 1
                  gxp(n) = rx * gxp(n-dj) + gxp(n-dj+di)
                  gyp(n) = ry * gyp(n-dj) + gyp(n-dj+di)
                  gzp(n) = rz * gzp(n-dj) + gzp(n-dj+di)
               end do
            end do
         end do
      end do
   end subroutine cint_g0_il2d_4d

   subroutine cint_g0_ik2d_4d(g, envs)
      real(dp), intent(inout), target :: g(0:)
      ! Three rows bound once, as the C's DEF_GXYZ does; see cint_g0_2e_2d.
      real(dp), pointer, contiguous :: gxp(:), gyp(:), gzp(:)
      type(cint_env_vars), intent(in) :: envs
      integer  :: lk, lj, ll, nmax, mmax, nroots, i, j, k, l, ptr, n
      integer  :: di, dk, dl, dj, gx, gy, gz
      real(dp) :: rx, ry, rz

      lk = envs%lk_ceil; lj = envs%lj_ceil; ll = envs%ll_ceil
      if (lj == 0 .and. ll == 0) return
      nmax = envs%li_ceil + envs%lj_ceil
      mmax = envs%lk_ceil + envs%ll_ceil
      nroots = envs%nrys_roots
      di = envs%g_stride_i; dk = envs%g_stride_k
      dl = envs%g_stride_l; dj = envs%g_stride_j
      gx = 0; gy = envs%g_size; gz = envs%g_size*2
      gxp(0:) => g(gx:); gyp(0:) => g(gy:); gzp(0:) => g(gz:)

      rx = envs%rkrl(0); ry = envs%rkrl(1); rz = envs%rkrl(2)
      do l = 1, ll
         do k = 0, mmax - l
            do i = 0, nmax
               ptr = l*dl + k*dk + i*di
               do n = ptr, ptr + nroots - 1
                  gxp(n) = rx * gxp(n-dl) + gxp(n-dl+dk)
                  gyp(n) = ry * gyp(n-dl) + gyp(n-dl+dk)
                  gzp(n) = rz * gzp(n-dl) + gzp(n-dl+dk)
               end do
            end do
         end do
      end do

      rx = envs%rirj(0); ry = envs%rirj(1); rz = envs%rirj(2)
      do j = 1, lj
         do l = 0, ll
            do k = 0, lk
               ptr = j*dj + l*dl + k*dk
               do n = ptr, ptr + dk - di*j - 1
                  gxp(n) = rx * gxp(n-dj) + gxp(n-dj+di)
                  gyp(n) = ry * gyp(n-dj) + gyp(n-dj+di)
                  gzp(n) = rz * gzp(n-dj) + gzp(n-dj+di)
               end do
            end do
         end do
      end do
   end subroutine cint_g0_ik2d_4d

   ! Rys roots for this primitive quartet, the b/c coefficients built from
   ! them, then the 2D recursion and the transfer into the 4D layout.
   ! Returns 0 when the quartet is negligible and nothing was written.
   function cint_g0_2e(g, rij, rkl, cutoff, envs) result(ok)
      real(dp), intent(inout) :: g(0:)
      ! EXPLICIT SHAPE, not assumed.  These are always the three components
      ! of a pair centre, and an assumed-shape dummy makes every read of
      ! them go through a descriptor: against the C, rij(0)-rkl(0) cost 5
      ! instructions to its 2, and rij(2)-rkl(2) cost 10 to its 3.
      real(dp), intent(in)    :: rij(0:2), rkl(0:2)
      real(dp), intent(in)    :: cutoff
      type(cint_env_vars), intent(in) :: envs
      integer :: ok

      integer  :: irys, nroots, rorder, err, w
      real(dp) :: aij, akl, a0, a1, fac1, x, omega, theta, sqrt_theta, ut
      real(dp) :: xij_kl, yij_kl, zij_kl, rr
      real(dp) :: u(0:MXRYSROOTS-1)
      real(dp) :: u2, tmp1, tmp2, tmp3, tmp4, tmp5
      real(dp) :: rijrx, rijry, rijrz, rklrx, rklry, rklrz
      type(rys_2e_t) :: bc

      ok = 0
      nroots = envs%nrys_roots
      w = envs%g_size * 2          ! the weights live in the gz block
      aij = envs%ai + envs%aj
      akl = envs%ak + envs%al
      xij_kl = rij(0) - rkl(0)
      yij_kl = rij(1) - rkl(1)
      zij_kl = rij(2) - rkl(2)
      rr = xij_kl*xij_kl + yij_kl*yij_kl + zij_kl*zij_kl

      a1 = aij * akl
      a0 = a1 / (aij + akl)
      fac1 = sqrt(a0 / (a1*a1*a1)) * envs%fac
      x = a0 * rr
      omega = envs%omega
      theta = 0.0_dp

      if (omega == 0.0_dp) then
         err = cint_rys_roots_lr(nroots, x, u, g(w:))
         if (err /= 0) return
      else if (omega < 0.0_dp) then
         ! short-range part of the range-separated Coulomb kernel
         theta = omega*omega / (omega*omega + a0)
         ! a vanishing erfc gives ~0 weights, which destabilises sr_rys_roots
         if (theta * x > cutoff .or. theta * x > EXPCUTOFF_SR) return
         rorder = envs%rys_order
         if (rorder == nroots) then
            err = cint_rys_roots_sr(nroots, x, sqrt(theta), u, g(w:))
            if (err /= 0) return
         else
            sqrt_theta = -sqrt(theta)
            err = cint_rys_roots_lr(rorder, x, u, g(w:))
            err = cint_rys_roots_lr(rorder, theta*x, u(rorder:), g(w+rorder:))
            if (envs%g_size == 2) then
               g(0) = 1.0_dp
               g(1) = 1.0_dp
               g(2) = 1.0_dp
               g(3) = 1.0_dp
               g(4) = g(4) * fac1
               g(5) = g(5) * fac1 * sqrt_theta
               ok = 1
               return
            end if
            do irys = rorder, nroots - 1
               ut = u(irys) * theta
               u(irys) = ut / (u(irys) + 1.0_dp - ut)
               g(w+irys) = g(w+irys) * sqrt_theta
            end do
         end if
      else
         theta = omega*omega / (omega*omega + a0)
         x = x * theta
         fac1 = fac1 * sqrt(theta)
         err = cint_rys_roots_lr(nroots, x, u, g(w:))
         if (err /= 0) return
         ! rescale u so the rest of the code is unchanged
         do irys = 0, nroots - 1
            ut = u(irys) * theta
            u(irys) = ut / (u(irys) + 1.0_dp - ut)
         end do
      end if

      if (envs%g_size == 1) then
         g(0) = 1.0_dp
         g(1) = 1.0_dp
         g(2) = g(2) * fac1
         ok = 1
         return
      end if

      rijrx = rij(0) - envs%rx_in_rijrx(0)
      rijry = rij(1) - envs%rx_in_rijrx(1)
      rijrz = rij(2) - envs%rx_in_rijrx(2)
      rklrx = rkl(0) - envs%rx_in_rklrx(0)
      rklry = rkl(1) - envs%rx_in_rklrx(1)
      rklrz = rkl(2) - envs%rx_in_rklrx(2)

      do irys = 0, nroots - 1
         u2 = a0 * u(irys)
         tmp4 = 0.5_dp / (u2 * (aij + akl) + a1)
         tmp5 = u2 * tmp4
         tmp1 = 2.0_dp * tmp5
         tmp2 = tmp1 * akl
         tmp3 = tmp1 * aij
         bc%b00(irys) = tmp5
         bc%b10(irys) = tmp5 + tmp4 * akl
         bc%b01(irys) = tmp5 + tmp4 * aij
         bc%c00x(irys) = rijrx - tmp2 * xij_kl
         bc%c00y(irys) = rijry - tmp2 * yij_kl
         bc%c00z(irys) = rijrz - tmp2 * zij_kl
         bc%c0px(irys) = rklrx + tmp3 * xij_kl
         bc%c0py(irys) = rklry + tmp3 * yij_kl
         bc%c0pz(irys) = rklrz + tmp3 * zij_kl
         g(w+irys) = g(w+irys) * fac1
      end do

      select case (envs%f_g0_2d4d)
      case (G2D4D_UNROLLED)
         call cint_g0_2e_2d4d_unrolled(g, bc, envs)
      case (G2D4D_SR_UNROLLED)
         call cint_srg0_2e_2d4d_unrolled(g, bc, envs)
      case (G2D4D_LJ)
         call cint_g0_2e_2d(g, bc, envs); call cint_g0_lj2d_4d(g, envs)
      case (G2D4D_KJ)
         call cint_g0_2e_2d(g, bc, envs); call cint_g0_kj2d_4d(g, envs)
      case (G2D4D_IL)
         call cint_g0_2e_2d(g, bc, envs); call cint_g0_il2d_4d(g, envs)
      case (G2D4D_IK)
         call cint_g0_2e_2d(g, bc, envs); call cint_g0_ik2d_4d(g, envs)
      case (G2D4D_NONE)
         call cint_g0_2e_2d(g, bc, envs)
      case default
         error stop "cint_g0_2e: no 2d4d transfer selected"
      end select
      ok = 1
   end function cint_g0_2e

   ! ---- the optimizer builder for this arity ----------------------------
   !
   ! It lives here rather than in cint_opt because it names this module's
   ! envs init, and cint_opt has to stay underneath every arity.

   subroutine cint_all_2e_optimizer(opt, ng, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: ng(0:), natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
      call opt_setij(opt, ng, atm, bas, nbas, env)
      call opt_set_non0coeff(opt, bas, nbas, env)
      call opt_gen_idx(opt, cint_init_int2e_envvars, cint_g2e_index_xyz, 4, 6, &
                       ng, atm, natm, bas, nbas, env)
      call opt_finish(opt, 4, ng, nbas)
   end subroutine cint_all_2e_optimizer

   subroutine cint_all_3c2e_optimizer(opt, ng, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: ng(0:), natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
      call opt_setij(opt, ng, atm, bas, nbas, env)
      call opt_set_non0coeff(opt, bas, nbas, env)
      call opt_gen_idx(opt, cint_init_int3c2e_envvars, cint_g2e_index_xyz, 3, 12, &
                       ng, atm, natm, bas, nbas, env)
      call opt_finish(opt, 3, ng, nbas)
   end subroutine cint_all_3c2e_optimizer

   subroutine cint_all_2c2e_optimizer(opt, ng, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: ng(0:), natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
      call opt_set_log_maxc(opt, bas, nbas, env)
      call opt_set_non0coeff(opt, bas, nbas, env)
      call opt_gen_idx(opt, cint_init_int2c2e_envvars, cint_g1e_index_xyz, 2, ANG_MAX, &
                       ng, atm, natm, bas, nbas, env)
      call opt_finish(opt, 2, ng, nbas)
   end subroutine cint_all_2c2e_optimizer


end module cint_g2e
