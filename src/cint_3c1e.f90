!
! Three-centre one-electron integrals: (ij|k) with no 1/r12 between them.
!
! Ported from src/g3c1e.c and src/cint3c1e.c.  Two operators share the file:
! the plain overlap <i j k>, and the nuclear-attraction form <i j| 1/r |k>,
! which is the same recursion evaluated at each Rys root of the 1/r kernel.
!
! THE G ARRAY IS BUILT ALONG j, NOT i.  Every other family in this port raises
! the bra index first; here the C runs the one-dimensional recursion up the j
! axis with `dj = li + 1`, transfers to i, and only then transfers to k -- and
! the j stride *changes* partway through, from li+1 to envs%g_stride_j.  That
! is not a simplification waiting to happen: the first phase works in a
! compressed (i,j) plane and the second in the full (i,j,k) block, and the two
! strides are the same number only when lk is zero.
!
module cint_3c1e
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws, ws_ensure, ws_alloc_d, ws_alloc_i, &
                             ws_mark, ws_rewind, &
                             ws_opt_log_maxc, ws_opt_non0, ws_opt_idx
   use cint_screen,    only: cint_non0coeff_byshell
   use cint_g1e,       only: cint_common_fac_sp, cint_prim_to_ctr_0, &
                             cint_prim_to_ctr_1, cint_nuc_mod
   use cint_rys_roots, only: cint_rys_roots_lr
   use cint_g2e,       only: cint_g2e_index_xyz
   use cint_bas,       only: cint_cart_comp, cint_square_dist, &
                             ANG_OF, NCTR_OF, ATOM_OF, NPRIM_OF, &
                             PTR_EXP, PTR_COEFF, PTR_COORD, CHARGE_OF, &
                             POINT_NUC, FRAC_CHARGE_NUC, PTR_FRAC_CHARGE, &
                             NUC_MOD_OF
   use cint_blas,      only: cint_dmat_transpose, cint_dplus_transpose
   use cint_3c2e,      only: apply_c2s_cart_3c2e1_pub => apply_c2s_cart_3c2e1, &
                             apply_c2s_sph_3c2e1_pub  => apply_c2s_sph_3c2e1, &
                             apply_c2s_dset0_3c_pub   => apply_c2s_dset0_3c
   use cint_opt,       only: cint_del_optimizer, opt_set_log_maxc, &
                             opt_set_non0coeff, opt_setij, opt_gen_idx, opt_finish
   implicit none
   private

   public :: cint_init_int3c1e_envvars, cint_g3c1e_index_xyz
   public :: cint_g3c1e_ovlp, cint_g3c1e_nuc, cint_gout3c1e
   public :: cint_3c1e_drv, int3c1e_cart, int3c1e_sph
   public :: int3c1e_rinv_cart, int3c1e_rinv_sph
   public :: int3c1e_optimizer
   public :: int3c1e_rinv_optimizer
   public :: cint_all_3c1e_optimizer

   ! which operator sits between the pair and the auxiliary shell
   integer, parameter, public :: INT3C1E_OVLP = 0
   integer, parameter, public :: INT3C1E_RINV = 1
   integer, parameter, public :: INT3C1E_NUC  = 2

   real(dp), parameter :: SQRTPI = 1.7724538509055160272981674833411451_dp
   real(dp), parameter :: PI     = 3.1415926535897932384626433832795029_dp

contains

   subroutine cint_init_int3c1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      type(cint_env_vars), intent(inout) :: envs
      integer,  intent(in) :: ng(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      integer :: i_sh, j_sh, k_sh, dli, dlj, dlk, nmax, ip

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
      ! note the order: i, j, k -- unlike the 2e families, which are (i,k,l,j)
      envs%nf = envs%nfi * envs%nfj * envs%nfk

      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, i_sh))
      envs%ri(0) = env(ip); envs%ri(1) = env(ip+1); envs%ri(2) = env(ip+2)
      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, j_sh))
      envs%rj(0) = env(ip); envs%rj(1) = env(ip+1); envs%rj(2) = env(ip+2)
      ip = atm_of(envs, PTR_COORD, bas_of(envs, ATOM_OF, k_sh))
      envs%rk(0) = env(ip); envs%rk(1) = env(ip+1); envs%rk(2) = env(ip+2)

      envs%gbits = ng(GSHIFT)
      envs%ncomp_e1 = ng(POS_E1)
      envs%ncomp_e2 = 0
      envs%ncomp_tensor = ng(TENSOR)

      envs%li_ceil = envs%i_l + ng(IINC)
      envs%lj_ceil = envs%j_l + ng(JINC)
      envs%lk_ceil = envs%k_l + ng(KINC)
      envs%ll_ceil = 0
      envs%nrys_roots = (envs%li_ceil + envs%lj_ceil + envs%lk_ceil)/2 + 1

      envs%common_factor = SQRTPI * PI &
         * cint_common_fac_sp(envs%i_l) * cint_common_fac_sp(envs%j_l) &
         * cint_common_fac_sp(envs%k_l)
      if (env(PTR_EXPCUTOFF) == 0.0_dp) then
         envs%expcutoff = EXPCUTOFF
      else
         envs%expcutoff = max(MIN_EXPCUTOFF, env(PTR_EXPCUTOFF))
      end if

      dli = envs%li_ceil + 1
      dlj = envs%lj_ceil + envs%lk_ceil + 1
      dlk = envs%lk_ceil + 1
      envs%g_stride_i = 1
      envs%g_stride_j = dli
      envs%g_stride_k = dli * dlj
      envs%g_stride_l = envs%g_stride_k
      nmax = envs%li_ceil + dlj
      ! the first phase needs dli*nmax, the second dli*dlj*dlk; whichever is
      ! larger has to fit
      envs%g_size = max(dli*dlj*dlk, dli*nmax)

      envs%rirj = envs%ri - envs%rj
   end subroutine cint_init_int3c1e_envvars

   subroutine cint_g3c1e_index_xyz(idx, envs)
      integer, intent(out) :: idx(0:)
      type(cint_env_vars), intent(in) :: envs
      integer :: i, j, k, n, dj, dk, ofx, ofy, ofz
      integer :: ofjx, ofjy, ofjz, ofkx, ofky, ofkz
      integer :: i_nx(0:CART_MAX-1), i_ny(0:CART_MAX-1), i_nz(0:CART_MAX-1)
      integer :: j_nx(0:CART_MAX-1), j_ny(0:CART_MAX-1), j_nz(0:CART_MAX-1)
      integer :: k_nx(0:CART_MAX-1), k_ny(0:CART_MAX-1), k_nz(0:CART_MAX-1)

      dj = envs%g_stride_j
      dk = envs%g_stride_k
      call cint_cart_comp(i_nx, i_ny, i_nz, envs%i_l)
      call cint_cart_comp(j_nx, j_ny, j_nz, envs%j_l)
      call cint_cart_comp(k_nx, k_ny, k_nz, envs%k_l)

      ofx = 0; ofy = envs%g_size; ofz = envs%g_size * 2
      n = 0
      do k = 0, envs%nfk - 1
         ofkx = ofx + dk * k_nx(k)
         ofky = ofy + dk * k_ny(k)
         ofkz = ofz + dk * k_nz(k)
         do j = 0, envs%nfj - 1
            ofjx = ofkx + dj * j_nx(j)
            ofjy = ofky + dj * j_ny(j)
            ofjz = ofkz + dj * j_nz(j)
            do i = 0, envs%nfi - 1
               idx(n+0) = ofjx + i_nx(i)
               idx(n+1) = ofjy + i_ny(i)
               idx(n+2) = ofjz + i_nz(i)
               n = n + 3
            end do
         end do
      end do
   end subroutine cint_g3c1e_index_xyz

   ! Kept because g3c1e.c defines it, and dropping it would quietly lose a
   ! piece of the C.  Nothing calls it there either: both loops in cint3c1e.c
   ! use CINTg2e_index_xyz instead.
   ! The plain overlap.  See the module note on the two j strides.
   subroutine cint_g3c1e_ovlp(g, og, ai, aj, ak, envs)
      ! TARGET plus one pointer per coordinate block, and the two
      ! three-vectors kept as scalars rather than array expressions -- see
      ! the note in cint_g1e_nuc.  This runs once per primitive triple.
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: og
      real(dp), intent(in)    :: ai, aj, ak
      type(cint_env_vars), intent(in) :: envs
      integer  :: li, lj, lk, nmax, mmax, gx, gy, gz, gs, dj, dk, i, j, k, off
      real(dp) :: aijk, aijk1, rjrk0, rjrk1, rjrk2, r0, r1, r2
      real(dp) :: rirj0, rirj1, rirj2
      real(dp), pointer, contiguous :: gxp(:), gyp(:), gzp(:)

      li = envs%li_ceil; lj = envs%lj_ceil; lk = envs%lk_ceil
      nmax = li + lj + lk
      mmax = lj + lk
      gs = envs%g_size
      gx = og; gy = og + gs; gz = og + gs*2
      gxp(0:) => g(gx:); gyp(0:) => g(gy:); gzp(0:) => g(gz:)
      gxp(0) = 1.0_dp
      gyp(0) = 1.0_dp
      gzp(0) = envs%fac
      if (nmax == 0) return

      dj = li + 1
      dk = envs%g_stride_k
      aijk = ai + aj + ak
      aijk1 = 0.5_dp / aijk
      rjrk0 = envs%rj(0) - envs%rk(0)
      rjrk1 = envs%rj(1) - envs%rk(1)
      rjrk2 = envs%rj(2) - envs%rk(2)
      ! division, not a reciprocal multiply: x/aijk and x*(1/aijk) round
      ! differently and the port has to match the C bit for bit
      r0 = envs%rj(0) - (ai*envs%ri(0) + aj*envs%rj(0) + ak*envs%rk(0)) / aijk
      r1 = envs%rj(1) - (ai*envs%ri(1) + aj*envs%rj(1) + ak*envs%rk(1)) / aijk
      r2 = envs%rj(2) - (ai*envs%ri(2) + aj*envs%rj(2) + ak*envs%rk(2)) / aijk
      rirj0 = envs%rirj(0); rirj1 = envs%rirj(1); rirj2 = envs%rirj(2)

      gxp(dj) = -r0 * gxp(0)
      gyp(dj) = -r1 * gyp(0)
      gzp(dj) = -r2 * gzp(0)

      do j = 1, nmax - 1
         gxp((j+1)*dj) = aijk1*j*gxp((j-1)*dj) - r0*gxp(j*dj)
         gyp((j+1)*dj) = aijk1*j*gyp((j-1)*dj) - r1*gyp(j*dj)
         gzp((j+1)*dj) = aijk1*j*gzp((j-1)*dj) - r2*gzp(j*dj)
      end do

      do i = 1, li
         do j = 0, nmax - i
            gxp(i+j*dj) = gxp(i-1+(j+1)*dj) - rirj0*gxp(i-1+j*dj)
            gyp(i+j*dj) = gyp(i-1+(j+1)*dj) - rirj1*gyp(i-1+j*dj)
            gzp(i+j*dj) = gzp(i-1+(j+1)*dj) - rirj2*gzp(i-1+j*dj)
         end do
      end do

      ! and now the full block, with the other j stride
      dj = envs%g_stride_j
      do k = 1, lk
         do j = 0, mmax - k
            off = k*dk + j*dj
            do i = off, off + li
               gxp(i) = gxp(i+dj-dk) + rjrk0*gxp(i-dk)
               gyp(i) = gyp(i+dj-dk) + rjrk1*gyp(i-dk)
               gzp(i) = gzp(i+dj-dk) + rjrk2*gzp(i-dk)
            end do
         end do
      end do
   end subroutine cint_g3c1e_ovlp

   ! The same recursion at one Rys root of the 1/r kernel.  t2 is that root
   ! transformed, and cr the charge centre.
   subroutine cint_g3c1e_nuc(g, og, ai, aj, ak, rijk, cr, t2, envs)
      ! Same treatment as cint_g3c1e_ovlp above.
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: og
      real(dp), intent(in)    :: ai, aj, ak, rijk(0:2), cr(0:2), t2
      type(cint_env_vars), intent(in) :: envs
      integer  :: li, lj, lk, nmax, mmax, gx, gy, gz, gs, dj, dk, i, j, k, off
      real(dp) :: aijk, aijk1, rjrk0, rjrk1, rjrk2, r0, r1, r2
      real(dp) :: rirj0, rirj1, rirj2
      real(dp), pointer, contiguous :: gxp(:), gyp(:), gzp(:)

      li = envs%li_ceil; lj = envs%lj_ceil; lk = envs%lk_ceil
      nmax = li + lj + lk
      mmax = lj + lk
      gs = envs%g_size
      gx = og; gy = og + gs; gz = og + gs*2
      gxp(0:) => g(gx:); gyp(0:) => g(gy:); gzp(0:) => g(gz:)
      gxp(0) = 1.0_dp
      gyp(0) = 1.0_dp
      gzp(0) = 2.0_dp/SQRTPI * envs%fac
      if (nmax == 0) return

      dj = li + 1
      dk = envs%g_stride_k
      aijk = ai + aj + ak
      rjrk0 = envs%rj(0) - envs%rk(0)
      rjrk1 = envs%rj(1) - envs%rk(1)
      rjrk2 = envs%rj(2) - envs%rk(2)
      r0 = envs%rj(0) - (rijk(0) + t2 * (cr(0) - rijk(0)))
      r1 = envs%rj(1) - (rijk(1) + t2 * (cr(1) - rijk(1)))
      r2 = envs%rj(2) - (rijk(2) + t2 * (cr(2) - rijk(2)))
      rirj0 = envs%rirj(0); rirj1 = envs%rirj(1); rirj2 = envs%rirj(2)

      gxp(dj) = -r0 * gxp(0)
      gyp(dj) = -r1 * gyp(0)
      gzp(dj) = -r2 * gzp(0)

      aijk1 = 0.5_dp * (1.0_dp - t2) / aijk
      do j = 1, nmax - 1
         gxp((j+1)*dj) = aijk1*j*gxp((j-1)*dj) - r0*gxp(j*dj)
         gyp((j+1)*dj) = aijk1*j*gyp((j-1)*dj) - r1*gyp(j*dj)
         gzp((j+1)*dj) = aijk1*j*gzp((j-1)*dj) - r2*gzp(j*dj)
      end do

      do i = 1, li
         do j = 0, nmax - i
            gxp(i+j*dj) = gxp(i-1+(j+1)*dj) - rirj0*gxp(i-1+j*dj)
            gyp(i+j*dj) = gyp(i-1+(j+1)*dj) - rirj1*gyp(i-1+j*dj)
            gzp(i+j*dj) = gzp(i-1+(j+1)*dj) - rirj2*gzp(i-1+j*dj)
         end do
      end do

      dj = envs%g_stride_j
      do k = 1, lk
         do j = 0, mmax - k
            off = k*dk + j*dj
            do i = off, off + li
               gxp(i) = gxp(i+dj-dk) + rjrk0*gxp(i-dk)
               gyp(i) = gyp(i+dj-dk) + rjrk1*gyp(i-dk)
               gzp(i) = gzp(i+dj-dk) + rjrk2*gzp(i-dk)
            end do
         end do
      end do
   end subroutine cint_g3c1e_nuc

   ! No Rys sum here: the 1/r root loop is outside, in the driver, because
   ! each root needs its own g array.
   subroutine cint_gout3c1e(gout, g, idx, envs, gout_empty)
      real(dp), intent(inout) :: gout(0:*)
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: idx(0:*)
      type(cint_env_vars), intent(in) :: envs
      integer,  intent(in)    :: gout_empty
      integer :: n, ix, iy, iz
      if (gout_empty /= 0) then
         do n = 0, envs%nf - 1
            ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
            gout(n) = g(ix) * g(iy) * g(iz)
         end do
      else
         do n = 0, envs%nf - 1
            ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
            gout(n) = gout(n) + g(ix) * g(iy) * g(iz)
         end do
      end if
   end subroutine cint_gout3c1e

   ! ---- the loop -------------------------------------------------------

   function cint_3c1e_loop(gctr, gctroff, envs, ws, fac, nuc_id, int_type, empty) &
         result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer,  intent(in)    :: gctroff
      type(cint_env_vars), intent(inout) :: envs
      ! TARGET for the block pointers bound below
      type(cint_ws), intent(inout), target :: ws
      real(dp), intent(in)    :: fac
      integer,  intent(in)    :: nuc_id, int_type
      logical,  intent(inout) :: empty
      logical :: has_value
      real(dp) :: facc

      integer :: i_sh, j_sh, k_sh, i_ctr, j_ctr, k_ctr, i_prim, j_prim, k_prim
      logical :: use_opt
      logical :: idx_cached
      integer :: ai_o, aj_o, ak_o, ci_o, cj_o, ck_o, n_comp, nf, nc
      integer :: leng, lenk, lenj, leni, len0
      integer :: og, ogout, ogctri, ogctrj, ogctrk, oidx, onon, osrt
      integer :: ip, jp, kp, i, err
      ! Hoisted as in the other loops; here the innermost body names the
      ! three exponents six times over, so it is the worst of them.
      real(dp), pointer, contiguous :: aip(:), ajp(:), akp(:)
      real(dp), pointer, contiguous :: cip(:), cjp(:), ckp(:)
      real(dp), pointer, contiguous :: gp(:), goutp(:)
      integer,  pointer, contiguous :: idxp(:), nonp(:)
      real(dp) :: aiv, ajv, akv
      real(dp) :: fac1i, fac1j, fac1k, eijk, dijk, aijk
      real(dp) :: aiajrr, aiakrr, ajakrr, rr_ij, rr_ik, rr_jk
      real(dp) :: rirk(0:2), rjrk(0:2), rijk(0:2), cr(0:2), tau, x, t2
      real(dp) :: u(0:MXRYSROOTS-1), w(0:MXRYSROOTS-1)
      integer, parameter :: E_I = 0, E_J = 1, E_K = 2, E_G = 3, E_M = 4
      logical :: flag(0:4)
      integer :: ii, ij, ik, ig, im, rys_empty

      has_value = .false.
      i_sh = envs%shls(0); j_sh = envs%shls(1); k_sh = envs%shls(2)
      i_ctr = envs%x_ctr(0); j_ctr = envs%x_ctr(1); k_ctr = envs%x_ctr(2)
      i_prim = bas_of(envs, NPRIM_OF, i_sh)
      j_prim = bas_of(envs, NPRIM_OF, j_sh)
      k_prim = bas_of(envs, NPRIM_OF, k_sh)
      ai_o = bas_of(envs, PTR_EXP, i_sh);   aj_o = bas_of(envs, PTR_EXP, j_sh)
      ak_o = bas_of(envs, PTR_EXP, k_sh)
      ci_o = bas_of(envs, PTR_COEFF, i_sh); cj_o = bas_of(envs, PTR_COEFF, j_sh)
      ck_o = bas_of(envs, PTR_COEFF, k_sh)
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor
      nf = envs%nf

      call ws_alloc_i(ws, nf*3, oidx)
      ! CINTg2e_index_xyz, not CINTg3c1e_index_xyz -- which g3c1e.c defines and
      ! cint3c1e.c never calls.  The two order the Cartesian components
      ! differently, (i,k,j) against (i,j,k), and they agree only when every
      ! shell is an s.  Following the C's actual call rather than its
      ! apparently-matching name is the whole difference.
      use_opt = .false.
      if (associated(ws%opt)) use_opt = cint_opt_usable(ws%opt, 3, envs)
      ! NOT `use_opt .and. ws_opt_idx(...)`: Fortran does not
      ! guarantee short-circuit evaluation, so the compiler is free to
      ! call ws_opt_idx even when use_opt is false -- and that
      ! dereferences ws%opt, which is null precisely when use_opt is
      ! false.  gfortran short-circuits at -O2 and does not at -O0, so
      ! this segfaulted only in a debug build.
      idx_cached = .false.
      if (use_opt) idx_cached = ws_opt_idx(ws, oidx, cint_opt_idx_key(envs, 3), nf*3)
      if (.not. idx_cached) then
         call cint_g2e_index_xyz(ws%i(oidx:), envs)
      end if
      call ws_alloc_i(ws, i_prim + j_prim + k_prim, onon)
      call ws_alloc_i(ws, i_prim*i_ctr + j_prim*j_ctr + k_prim*k_ctr, osrt)
      if (use_opt) then
         call ws_opt_non0(ws, onon, osrt, i_sh, i_prim, i_ctr, 0, 0)
         call ws_opt_non0(ws, onon, osrt, j_sh, j_prim, j_ctr, i_prim, i_prim*i_ctr)
         call ws_opt_non0(ws, onon, osrt, k_sh, k_prim, k_ctr, i_prim+j_prim, &
                          i_prim*i_ctr+j_prim*j_ctr)
      else
      call cint_non0coeff_byshell(ws%i(osrt:), ws%i(onon:), envs%env, ci_o, i_prim, i_ctr)
      call cint_non0coeff_byshell(ws%i(osrt+i_prim*i_ctr:), ws%i(onon+i_prim:), &
                                  envs%env, cj_o, j_prim, j_ctr)
      call cint_non0coeff_byshell(ws%i(osrt+i_prim*i_ctr+j_prim*j_ctr:), &
                                  ws%i(onon+i_prim+j_prim:), envs%env, ck_o, k_prim, k_ctr)
      end if

      nc = i_ctr * j_ctr * k_ctr
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenk = nf * nc * n_comp
      lenj = nf * i_ctr * j_ctr * n_comp
      leni = nf * i_ctr * n_comp
      len0 = nf * n_comp

      call ws_alloc_d(ws, leng, og)
      flag = .true.
      flag(E_M) = empty
      ii = E_I; ij = E_J; ik = E_K; ig = E_G; im = E_M
      if (n_comp == 1) then
         ogctrk = gctroff; ik = im
      else
         call ws_alloc_d(ws, lenk, ogctrk)
      end if
      if (k_ctr == 1) then
         ogctrj = ogctrk; ij = ik
      else
         call ws_alloc_d(ws, lenj, ogctrj)
      end if
      if (j_ctr == 1) then
         ogctri = ogctrj; ii = ij
      else
         call ws_alloc_d(ws, leni, ogctri)
      end if
      if (i_ctr == 1) then
         ogout = ogctri; ig = ii
      else
         call ws_alloc_d(ws, len0, ogout)
      end if

      ! The nuc loop folds the common factor into its own `fac` argument --
      ! one line, `fac *= envs->common_factor`, before the primitive loops.
      ! Missing it makes the answer wrong by exactly that factor, which for
      ! three s shells is a clean 8.
      facc = fac * envs%common_factor

      aip(0:) => envs%env(ai_o:); ajp(0:) => envs%env(aj_o:)
      akp(0:) => envs%env(ak_o:)
      cip(0:) => envs%env(ci_o:); cjp(0:) => envs%env(cj_o:)
      ckp(0:) => envs%env(ck_o:)
      gp(0:)    => ws%d(og:og+leng-1)
      goutp(0:) => ws%d(ogout:ogout+len0-1)
      idxp(0:)  => ws%i(oidx:oidx+nf*3-1)
      nonp(0:)  => ws%i(onon:onon+i_prim+j_prim+k_prim-1)

      rirk = envs%ri - envs%rk
      rjrk = envs%rj - envs%rk
      rr_ij = sum(envs%rirj * envs%rirj)
      rr_ik = sum(rirk * rirk)
      rr_jk = sum(rjrk * rjrk)

      if (int_type /= INT3C1E_OVLP) then
         if (nuc_id < 0) then
            cr = envs%env(PTR_RINV_ORIG:PTR_RINV_ORIG+2)
         else
            i = atm_of(envs, PTR_COORD, nuc_id)
            cr = envs%env(i:i+2)
         end if
      end if

      do kp = 0, k_prim - 1
         akv = akp(kp)
         envs%ak = akv
         if (k_ctr == 1) then
            if (int_type == INT3C1E_OVLP) then
               fac1k = envs%common_factor * ckp(kp)
            else
               fac1k = facc * ckp(kp)
            end if
         else
            if (int_type == INT3C1E_OVLP) then
               fac1k = envs%common_factor
            else
               fac1k = facc
            end if
            flag(ij) = .true.
         end if

         do jp = 0, j_prim - 1
            ajv = ajp(jp)
            envs%aj = ajv
            if (j_ctr == 1) then
               fac1j = fac1k * cjp(jp)
            else
               fac1j = fac1k
               flag(ii) = .true.
            end if
            ajakrr = ajv * akv * rr_jk
            do ip = 0, i_prim - 1
               aiv = aip(ip)
               envs%ai = aiv
               aijk = aiv + ajv + akv
               aiakrr = aiv * akv * rr_ik
               aiajrr = aiv * ajv * rr_ij
               eijk = (aiajrr + aiakrr + ajakrr) / aijk
               if (eijk > EXPCUTOFF) cycle
               if (i_ctr == 1) then
                  fac1i = fac1j * cip(ip) * exp(-eijk)
               else
                  fac1i = fac1j * exp(-eijk)
               end if

               if (int_type == INT3C1E_OVLP) then
                  dijk = fac1i / (aijk * sqrt(aijk))
                  envs%fac = dijk
                  call cint_g3c1e_ovlp(ws%d, og, aiv, ajv, akv, envs)
                  call envs%f_gout(goutp, gp, idxp, envs, merge(1, 0, flag(ig)))
               else
                  dijk = fac1i / aijk
                  rijk = (aiv*envs%ri + ajv*envs%rj + akv*envs%rk) / aijk
                  tau = cint_nuc_mod(aijk, nuc_id, envs%atm, envs%env)
                  x = aijk * cint_square_dist(rijk, cr) * tau * tau
                  err = cint_rys_roots_lr(envs%nrys_roots, x, u, w)
                  rys_empty = merge(1, 0, flag(ig))
                  do i = 0, envs%nrys_roots - 1
                     t2 = u(i) / (1.0_dp + u(i)) * tau * tau
                     envs%fac = dijk * w(i) * tau
                     call cint_g3c1e_nuc(ws%d, og, aiv, ajv, akv, rijk, cr, t2, envs)
                     call envs%f_gout(goutp, gp, idxp, envs, rys_empty)
                     rys_empty = 0
                  end do
               end if
               if (i_ctr > 1) then
                  if (flag(ii)) then
                     call cint_prim_to_ctr_0(ws%d, ogctri, ogout, envs%env, &
                                             ci_o+ip, len0, i_prim, i_ctr)
                  else
                     call cint_prim_to_ctr_1(ws%d, ogctri, ogout, envs%env, &
                                             ci_o+ip, len0, i_prim, i_ctr, &
                                             nonp(ip), ws%i, osrt + ip*i_ctr)
                  end if
               end if
               flag(ii) = .false.
            end do
            if (.not. flag(ii)) then
               if (j_ctr > 1) then
                  if (flag(ij)) then
                     call cint_prim_to_ctr_0(ws%d, ogctrj, ogctri, envs%env, &
                                             cj_o+jp, leni, j_prim, j_ctr)
                  else
                     call cint_prim_to_ctr_1(ws%d, ogctrj, ogctri, envs%env, &
                                             cj_o+jp, leni, j_prim, j_ctr, &
                                             nonp(i_prim+jp), ws%i, &
                                             osrt + i_prim*i_ctr + jp*j_ctr)
                  end if
               end if
               flag(ij) = .false.
            end if
         end do
         if (.not. flag(ij)) then
            if (k_ctr > 1) then
               if (flag(ik)) then
                  call cint_prim_to_ctr_0(ws%d, ogctrk, ogctrj, envs%env, &
                                          ck_o+kp, lenj, k_prim, k_ctr)
               else
                  call cint_prim_to_ctr_1(ws%d, ogctrk, ogctrj, envs%env, &
                                          ck_o+kp, lenj, k_prim, k_ctr, &
                                          nonp(i_prim+j_prim+kp), ws%i, &
                                          osrt + i_prim*i_ctr + j_prim*j_ctr + kp*k_ctr)
               end if
            end if
            flag(ik) = .false.
         end if
      end do

      if (n_comp > 1 .and. .not. flag(ik)) then
         if (flag(im)) then
            call cint_dmat_transpose(gctr(gctroff:), ws%d(ogctrk:), nf*nc, n_comp)
            flag(im) = .false.
         else
            call cint_dplus_transpose(gctr(gctroff:), ws%d(ogctrk:), nf*nc, n_comp)
         end if
      end if
      empty = flag(im)
      has_value = .not. empty
   end function cint_3c1e_loop


   ! ---- driver ---------------------------------------------------------

   function cint_3c1e_drv(out, dims, envs, ws, c2s_kind, int_type) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      integer,  intent(in)    :: c2s_kind, int_type
      logical :: has_value

      integer :: nc, n_comp, ogctr, obuf, n, nout, counts(0:3), nd, ni
      integer :: ip, jp, kp, leng, len0, ia, dmark, imark
      real(dp) :: fac
      logical :: empty, dummy

      nc = envs%nf * envs%x_ctr(0) * envs%x_ctr(1) * envs%x_ctr(2)
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor

      ip = bas_of(envs, NPRIM_OF, envs%shls(0))
      jp = bas_of(envs, NPRIM_OF, envs%shls(1))
      kp = bas_of(envs, NPRIM_OF, envs%shls(2))
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      len0 = envs%nf * n_comp
      nd = leng + len0 + nc*n_comp*4 &
         + 4 * envs%nfi * envs%nfk * (2*envs%j_l + 1) + 128
      ni = envs%nf*3 + ip*envs%x_ctr(0) + jp*envs%x_ctr(1) + kp*envs%x_ctr(2) &
         + ip + jp + kp + 64
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*n_comp, ogctr)
      ! NOT zeroed.  No C driver zeroes gctr -- CINT1e_drv, CINT3c2e_drv,
      ! CINT2c2e_drv, CINT3c1e_drv and CINT1e_grids_drv all just take the
      ! buffer and thread an empty flag through the loop, so the first
      ! contraction sets rather than accumulates.  The port threads the same
      ! flag, so the memset was pure extra writes.

      empty = .true.
      select case (int_type)
      case (INT3C1E_OVLP)
         dummy = cint_3c1e_loop(ws%d, ogctr, envs, ws, 1.0_dp, -1, INT3C1E_OVLP, empty)
      case (INT3C1E_RINV)
         dummy = cint_3c1e_loop(ws%d, ogctr, envs, ws, 1.0_dp, -1, INT3C1E_RINV, empty)
      case default
         ! sum over nuclei, each with its own charge
         call ws_alloc_d(ws, nc*n_comp, obuf)
         do ia = 0, envs%natm - 1
            fac = -real(abs(atm_of(envs, CHARGE_OF, ia)), dp)
            if (atm_of(envs, NUC_MOD_OF, ia) == FRAC_CHARGE_NUC) then
               fac = -abs(envs%env(atm_of(envs, PTR_FRAC_CHARGE, ia)))
            end if
            dummy = cint_3c1e_loop(ws%d, ogctr, envs, ws, fac, ia, INT3C1E_NUC, empty)
         end do
      end select
      has_value = .not. empty

      select case (c2s_kind)
      case (C2S_SPH_3C2E1)
         counts(0) = (envs%i_l*2+1) * envs%x_ctr(0)
         counts(1) = (envs%j_l*2+1) * envs%x_ctr(1)
         counts(2) = (envs%k_l*2+1) * envs%x_ctr(2)
      case default
         counts(0) = envs%nfi * envs%x_ctr(0)
         counts(1) = envs%nfj * envs%x_ctr(1)
         counts(2) = envs%nfk * envs%x_ctr(2)
      end select
      counts(3) = 1
      nout = dims(0) * dims(1) * dims(2)

      if (has_value) then
         call ws_mark(ws, dmark, imark)
         do n = 0, n_comp - 1
            call ws_rewind(ws, dmark, imark)
            select case (c2s_kind)
            case (C2S_SPH_3C2E1)
               call apply_c2s_sph_3c2e1_pub(out(nout*n:), ogctr + nc*n, dims, envs, ws)
            case default
               call apply_c2s_cart_3c2e1_pub(out(nout*n:), ws%d(ogctr+nc*n:), dims, envs)
            end select
         end do
      else
         do n = 0, n_comp - 1
            call apply_c2s_dset0_3c_pub(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_3c1e_drv

   ! ---- entry points ---------------------------------------------------

   function int3c1e_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int3c1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout3c1e
      hv = cint_3c1e_drv(out, dims, envs, ws, C2S_CART_3C2E1, INT3C1E_OVLP)
   end function int3c1e_cart

   function int3c1e_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int3c1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout3c1e
      hv = cint_3c1e_drv(out, dims, envs, ws, C2S_SPH_3C2E1, INT3C1E_OVLP)
   end function int3c1e_sph

   function int3c1e_rinv_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int3c1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout3c1e
      hv = cint_3c1e_drv(out, dims, envs, ws, C2S_CART_3C2E1, INT3C1E_RINV)
   end function int3c1e_rinv_cart

   function int3c1e_rinv_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int3c1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout3c1e
      hv = cint_3c1e_drv(out, dims, envs, ws, C2S_SPH_3C2E1, INT3C1E_RINV)
   end function int3c1e_rinv_sph

   ! ---- optimizers ------------------------------------------------------

   subroutine int3c1e_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_all_3c1e_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int3c1e_optimizer

   subroutine int3c1e_rinv_optimizer(opt, atm, natm, bas, nbas, env)
      ! The C's is `*opt = NULL`: this family has no optimizer upstream.
      ! Clearing rather than ignoring, so a stale optimizer left on the
      ! workspace cannot be picked up by the next call.
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
   end subroutine int3c1e_rinv_optimizer


   ! ---- the optimizer builder for this arity ----------------------------
   !
   ! It lives here rather than in cint_opt because it names this module's
   ! envs init, and cint_opt has to stay underneath every arity.

   subroutine cint_all_3c1e_optimizer(opt, ng, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: ng(0:), natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      call cint_del_optimizer(opt)
      call opt_setij(opt, ng, atm, bas, nbas, env)
      call opt_set_non0coeff(opt, bas, nbas, env)
      call opt_gen_idx(opt, cint_init_int3c1e_envvars, cint_g2e_index_xyz, 3, 12, &
                       ng, atm, natm, bas, nbas, env)
      call opt_finish(opt, 3, ng, nbas)
   end subroutine cint_all_3c1e_optimizer


end module cint_3c1e
