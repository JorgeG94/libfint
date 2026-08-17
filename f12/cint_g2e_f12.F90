! Guarded by COMPILE_F12.  Without it this file compiles to nothing at all,
! which is what lets every source in fortran/src be listed unconditionally --
! and the reason it has to be a guard rather than "just always build it" is
! that the chain below ends at cint_tab_stg_roots, a 100 MB generated table
! that is deliberately not in the repository.  See
! scripts/stg_roots_to_fortran.py.
!
! No #else stub: fpm scans for module names lexically, before the
! preprocessor runs, so a second `module cint_g2e_f12` in this file reads to it
! as a duplicate definition.  An absent module is also the better error for
! anyone who USEs it without asking for F12.
#ifdef COMPILE_F12
!
! The F12 kernels: the Yukawa potential and the Slater-type geminal.
!
! Ported from src/g2e_f12.c and src/cint2e_f12.c.  Only built with WITH_F12,
! which is off by default here because it is off by default upstream.
!
! Both kernels are the ordinary two-electron machinery with two changes.  The
! prefactor is fac/(sqrt(aij+akl)*aij*akl) rather than sqrt(a0/a1^3)*fac, and
! when the geminal exponent zeta is positive the quadrature comes from
! CINTstg_roots rather than the Rys roots -- after which the roots and
! weights are rescaled, differently for the two kernels.  With zeta at or
! below zero both fall back to the Rys roots and the rescaling is skipped,
! which is how the C spells "no geminal".
!
! The envs setup differs from the ordinary one in three places, all of them
! the C's: no +1 on the exponent cutoff, ceil(L/2)+1 roots with no
! short-range doubling, and no unrolled 2D-to-4D arm.
!
module cint_g2e_f12
   use cint_const,     only: dp
   use cint_envs
   use cint_bas,       only: ANG_OF, NCTR_OF, ATOM_OF, PTR_COORD
   use cint_g1e,       only: cint_common_fac_sp
   use cint_g2e
   use cint_rys_roots, only: cint_rys_roots_lr
   use cint_stg_quad,  only: cint_stg_roots
   use cint_2e,        only: cint_2e_set_f12_hook
   use cint_opt,       only: cint_del_optimizer, opt_set_non0coeff, opt_setij, &
                             opt_gen_idx, opt_finish
   implicit none
   private

   real(dp), parameter :: PI     = 3.1415926535897932384626433832795029_dp
   real(dp), parameter :: SQRTPI = 1.7724538509055160272981674833411451_dp

   public :: cint_init_int2e_yp_envvars, cint_init_int2e_stg_envvars
   public :: cint_g0_2e_yp, cint_g0_2e_stg
   public :: cint_all_2e_stg_optimizer

   ! Which rescaling the geminal weights get.  The C says this with two
   ! nearly identical functions; one flag on the envs says it once.
   integer, parameter :: F12_YP = 0, F12_STG = 1

   ! Which of the two rescalings the registered hook applies.  A module
   ! variable rather than a field on the envs, for the same reason the hook
   ! is: the envs is declared on the stack of every entry point in the
   ! library, and this family is one of ten.
   logical, save :: f12_use_stg = .false.
   !$omp threadprivate(f12_use_stg)

contains

   subroutine cint_init_int2e_f12_envvars(envs, ng, shls, atm, natm, bas, nbas, env, kind)
      type(cint_env_vars), intent(inout) :: envs
      integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
      ! no INTENT(IN) on the three tables: envs points at them, and F2018
      ! 8.5.10 forbids an INTENT(IN) dummy as a pointer target.  They are
      ! never written -- see the note on cint_env_vars.
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      integer,  intent(in) :: kind

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

      i_sh = shls(0); j_sh = shls(1); k_sh = shls(2); l_sh = shls(3)
      envs%i_l = bas_of(envs, ANG_OF, i_sh)
      envs%j_l = bas_of(envs, ANG_OF, j_sh)
      envs%k_l = bas_of(envs, ANG_OF, k_sh)
      envs%l_l = bas_of(envs, ANG_OF, l_sh)
      envs%x_ctr(0) = bas_of(envs, NCTR_OF, i_sh)
      envs%x_ctr(1) = bas_of(envs, NCTR_OF, j_sh)
      envs%x_ctr(2) = bas_of(envs, NCTR_OF, k_sh)
      envs%x_ctr(3) = bas_of(envs, NCTR_OF, l_sh)
      envs%nfi = (envs%i_l+1)*(envs%i_l+2)/2
      envs%nfj = (envs%j_l+1)*(envs%j_l+2)/2
      envs%nfk = (envs%k_l+1)*(envs%k_l+2)/2
      envs%nfl = (envs%l_l+1)*(envs%l_l+2)/2
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
         * cint_common_fac_sp(envs%i_l) * cint_common_fac_sp(envs%j_l) &
         * cint_common_fac_sp(envs%k_l) * cint_common_fac_sp(envs%l_l)
      if (env(PTR_EXPCUTOFF) == 0.0_dp) then
         envs%expcutoff = EXPCUTOFF
      else
         ! no +1 here, unlike the ordinary two-electron setup: g2e_f12.c
         ! does not add it, and the screening has to match the C's
         envs%expcutoff = max(MIN_EXPCUTOFF, env(PTR_EXPCUTOFF))
      end if

      envs%gbits = ng(GSHIFT)
      envs%ncomp_e1 = ng(POS_E1)
      envs%ncomp_e2 = ng(POS_E2)
      envs%ncomp_tensor = ng(TENSOR)

      envs%li_ceil = envs%i_l + ng(IINC)
      envs%lj_ceil = envs%j_l + ng(JINC)
      envs%lk_ceil = envs%k_l + ng(KINC)
      envs%ll_ceil = envs%l_l + ng(LINC)
      ! ceil(L_tot/2) + 1, and no short-range doubling: the attenuation here
      ! is the geminal, not an erfc kernel.
      nrys_roots = (envs%li_ceil + envs%lj_ceil + envs%lk_ceil + envs%ll_ceil + 3)/2
      rys_order = nrys_roots
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

      ! No unrolled arm.  g2e_f12.c picks only the four general transfers,
      ! and its lj case names CINTg0_2e_stg_lj2d4d -- a byte-for-byte
      ! duplicate of CINTg0_2e_lj2d4d that exists only because g2e_f12.c
      ! cannot see the static one.  Same transfer, so G2D4D_LJ.
      if (kbase) then
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
      envs%g0_kind = G0_F12
      ! The kernel itself goes through the hook cint_2e holds, registered
      ! once below; the yp/stg choice rides on f12_use_stg because the two
      ! differ only in four lines of weight rescaling.
      f12_use_stg = (kind == F12_STG)
      call cint_2e_set_f12_hook(g0_f12_hook)
   end subroutine cint_init_int2e_f12_envvars


   subroutine cint_init_int2e_yp_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      type(cint_env_vars), intent(inout) :: envs
      integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      call cint_init_int2e_f12_envvars(envs, ng, shls, atm, natm, bas, nbas, env, F12_YP)
   end subroutine cint_init_int2e_yp_envvars

   subroutine cint_init_int2e_stg_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      type(cint_env_vars), intent(inout) :: envs
      integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      call cint_init_int2e_f12_envvars(envs, ng, shls, atm, natm, bas, nbas, env, F12_STG)
   end subroutine cint_init_int2e_stg_envvars

   function g0_f12_hook(g, rij, rkl, cutoff, envs) result(ok)
      real(dp), intent(inout) :: g(0:)
      real(dp), intent(in)    :: rij(0:2), rkl(0:2)
      real(dp), intent(in)    :: cutoff
      type(cint_env_vars), intent(in) :: envs
      integer :: ok
      if (f12_use_stg) then
         ok = f12_kernel(g, rij, rkl, cutoff, envs, F12_STG)
      else
         ok = f12_kernel(g, rij, rkl, cutoff, envs, F12_YP)
      end if
   end function g0_f12_hook

   function cint_g0_2e_yp(g, rij, rkl, cutoff, envs) result(ok)
      real(dp), intent(inout) :: g(0:)
      real(dp), intent(in)    :: rij(0:2), rkl(0:2)
      real(dp), intent(in)    :: cutoff
      type(cint_env_vars), intent(in) :: envs
      integer :: ok
      ok = f12_kernel(g, rij, rkl, cutoff, envs, F12_YP)
   end function cint_g0_2e_yp

   function cint_g0_2e_stg(g, rij, rkl, cutoff, envs) result(ok)
      real(dp), intent(inout) :: g(0:)
      real(dp), intent(in)    :: rij(0:2), rkl(0:2)
      real(dp), intent(in)    :: cutoff
      type(cint_env_vars), intent(in) :: envs
      integer :: ok
      ok = f12_kernel(g, rij, rkl, cutoff, envs, F12_STG)
   end function cint_g0_2e_stg

   ! CINTg0_2e_yp and CINTg0_2e_stg (src/g2e_f12.c) are the same function
   ! twice over, differing in four lines of weight rescaling.  One body here,
   ! with the difference where it actually is.
   function f12_kernel(g, rij, rkl, cutoff, envs, kind) result(ok)
      real(dp), intent(inout) :: g(0:)
      real(dp), intent(in)    :: rij(0:2), rkl(0:2)
      real(dp), intent(in)    :: cutoff
      type(cint_env_vars), intent(in) :: envs
      integer,  intent(in)    :: kind
      integer :: ok

      integer  :: irys, nroots, w, err
      real(dp) :: aij, akl, a0, a1, fac1, x, ua, ua2, zeta
      real(dp) :: xij_kl, yij_kl, zij_kl, rr, u2, tmp1, tmp2, tmp3, tmp4, tmp5
      real(dp) :: rijrx, rijry, rijrz, rklrx, rklry, rklrz
      real(dp) :: u(0:MXRYSROOTS-1)
      type(rys_2e_t) :: bc

      ok = 0
      w = envs%g_size * 2               ! the weights live in the gz block
      zeta = envs%env(PTR_F12_ZETA)
      nroots = envs%nrys_roots
      aij = envs%ai + envs%aj
      akl = envs%ak + envs%al
      a1 = aij * akl
      a0 = a1 / (aij + akl)
      ! not sqrt(a0/(a1*a1*a1))*fac, as the ordinary kernel has it
      fac1 = envs%fac / (sqrt(aij + akl) * a1)

      ua = 0.0_dp
      if (zeta > 0.0_dp) ua = 0.25_dp * zeta * zeta / a0

      xij_kl = rij(0) - rkl(0)
      yij_kl = rij(1) - rkl(1)
      zij_kl = rij(2) - rkl(2)
      rr = xij_kl*xij_kl + yij_kl*yij_kl + zij_kl*zij_kl
      x = a0 * rr

      if (zeta > 0.0_dp) then
         call cint_stg_roots(nroots, x, ua, u, g(w:))
      else
         ! The C ignores this return, and a failed solve then leaves garbage
         ! rather than a zero block.  Inherited: changing it would make the
         ! two disagree exactly where the C is already wrong.
         err = cint_rys_roots_lr(nroots, x, u, g(w:))
      end if

      if (zeta > 0.0_dp) then
         if (kind == F12_STG) then
            ua2 = 2.0_dp * ua / zeta
            do irys = 0, nroots - 1
               ! w * ((1-u)*ua2), not (w*(1-u))*ua2: the C multiplies the
               ! pair first, and the two associations differ in the last bit.
               g(w+irys) = g(w+irys) * ((1.0_dp - u(irys)) * ua2)
               u(irys) = u(irys) / (1.0_dp - u(irys))
            end do
         else
            do irys = 0, nroots - 1
               g(w+irys) = g(w+irys) * u(irys)
               u(irys) = u(irys) / (1.0_dp - u(irys))
            end do
         end if
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
      case (G2D4D_LJ)
         call cint_g0_2e_2d(g, bc, envs); call cint_g0_lj2d_4d(g, envs)
      case (G2D4D_KJ)
         call cint_g0_2e_2d(g, bc, envs); call cint_g0_kj2d_4d(g, envs)
      case (G2D4D_IL)
         call cint_g0_2e_2d(g, bc, envs); call cint_g0_il2d_4d(g, envs)
      case (G2D4D_IK)
         call cint_g0_2e_2d(g, bc, envs); call cint_g0_ik2d_4d(g, envs)
      case default
         error stop "f12_kernel: no 2d4d transfer selected"
      end select
      ok = 1
   end function f12_kernel

   ! CINTall_2e_stg_optimizer (src/optimizer.c, behind WITH_F12).  Same
   ! caches as the ordinary two-electron one; only the envs init differs,
   ! which is exactly why opt_gen_idx takes it as an argument.
   subroutine cint_all_2e_stg_optimizer(opt, ng, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: ng(0:), natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
      call opt_setij(opt, ng, atm, bas, nbas, env)
      call opt_set_non0coeff(opt, bas, nbas, env)
      call opt_gen_idx(opt, cint_init_int2e_stg_envvars, cint_g2e_index_xyz, &
                       4, 6, ng, atm, natm, bas, nbas, env)
      call opt_finish(opt, 4, ng, nbas)
   end subroutine cint_all_2e_stg_optimizer

end module cint_g2e_f12
#endif
