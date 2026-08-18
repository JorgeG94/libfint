!
! The two-electron driver, and int2e on top of it.
!
! Ported from src/cint2e.c.  Four nested primitive loops (l, k, j, i), each
! contracting into the next, with screening applied at three levels: the kl
! pair, the ij pair, and a per-quartet cutoff handed down to cint_g0_2e so the
! Rys evaluation itself can bail out.
!
! Only the unoptimised loop is here.  The C has five more specialisations of
! it, picked by which contraction counts are 1, plus a path that consults a
! prebuilt CINTOpt; those are D8, and this one produces the same numbers more
! slowly.
!
! The `empty` flags are aliased the way the C aliases them -- see the note in
! cint_1e.f90, where getting that wrong by reasoning it out cost an afternoon.
!
module cint_2e
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws, ws_ensure, ws_ensure_pd, ws_alloc_d, ws_alloc_i, &
                             ws_mark, ws_rewind, ws_opt_log_maxc, ws_opt_non0
   use cint_screen,    only: pair_data, cint_set_pairdata, &
                             cint_log_max_pgto_coeff, cint_non0coeff_byshell
   use cint_g1e,       only: cint_prim_to_ctr_0, cint_prim_to_ctr_1, &
                             cint_prim_to_ctr_sp_0, cint_prim_to_ctr_sp_1
   use cint_g2e
   use cint_bas,       only: NPRIM_OF, PTR_EXP, PTR_COEFF, NF_SP
   use cint_blas,      only: cint_dmat_transpose, cint_dplus_transpose
   use cint_cart2sph,  only: cint_c2s_bra_sph, cint_c2s_ket_sph, &
                                RESULT_IN_GCART, RESULT_IN_GSPH
   use cint_g2e,       only: cint_all_2e_optimizer
   implicit none
   private

   public :: cint_gout2e, cint_2e_loop_nopt, cint_2e_drv, int2e_cache_size
   public :: int2e_cart, int2e_sph
   public :: cint_2e_set_f12_hook

   ! The F12 kernels live behind WITH_F12, in a module this one cannot name
   ! without making the dependency conditional.  They register themselves
   ! here instead, which is the same indirection the C gets from a function
   ! pointer on the envs -- kept off the envs because a pointer component
   ! there costs 1.9% across all 620 entry points.
   procedure(g0_2e_iface), pointer, save :: f12_hook => null()
   public :: int2e_optimizer

   real(dp), parameter :: R_GUESS = 8.0_dp

contains

   subroutine cint_2e_set_f12_hook(f)
      procedure(g0_2e_iface) :: f
      f12_hook => f
   end subroutine cint_2e_set_f12_hook

   ! Contract the g intermediates over the Rys roots.
   !
   ! UNROLLED ON THE ROOT COUNT, as the C is.  This was a plain loop, with a
   ! comment saying the two were the same sum in the same order -- which is
   ! true, and is why the numbers did not change when this was rewritten.
   ! What the comment missed is that a loop whose trip count the compiler
   ! cannot see is a loop it cannot unroll, and VTune put 0.07s of a 0.23s
   ! function on the loop header alone.  Most quartets have one Rys root, and
   ! for those the C emits three multiplies and no loop at all.
   !
   ! The sums stay left to right, so the unrolled form is bit-identical to the
   ! accumulator it replaces: adding the first term to a fresh zero is exact.
   subroutine cint_gout2e(gout, g, idx, envs, gout_empty)
      real(dp), intent(inout) :: gout(0:*)
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: idx(0:*)
      type(cint_env_vars), intent(in) :: envs
      integer,  intent(in)    :: gout_empty
      integer  :: n, i, ix, iy, iz, nf, nroots
      real(dp) :: s

      nf = envs%nf
      nroots = envs%nrys_roots
      if (gout_empty /= 0) then
         select case (nroots)
         case (1)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = g(ix) * g(iy) * g(iz)
            end do
         case (2)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = g(ix)   * g(iy)   * g(iz) &
                       + g(ix+1) * g(iy+1) * g(iz+1)
            end do
         case (3)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = g(ix)   * g(iy)   * g(iz) &
                       + g(ix+1) * g(iy+1) * g(iz+1) &
                       + g(ix+2) * g(iy+2) * g(iz+2)
            end do
         case (4)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = g(ix)   * g(iy)   * g(iz) &
                       + g(ix+1) * g(iy+1) * g(iz+1) &
                       + g(ix+2) * g(iy+2) * g(iz+2) &
                       + g(ix+3) * g(iy+3) * g(iz+3)
            end do
         case (5)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = g(ix) * g(iy) * g(iz) &
                       + g(ix+1) * g(iy+1) * g(iz+1) &
                       + g(ix+2) * g(iy+2) * g(iz+2) &
                       + g(ix+3) * g(iy+3) * g(iz+3) &
                       + g(ix+4) * g(iy+4) * g(iz+4)
            end do
         case (6)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = g(ix) * g(iy) * g(iz) &
                       + g(ix+1) * g(iy+1) * g(iz+1) &
                       + g(ix+2) * g(iy+2) * g(iz+2) &
                       + g(ix+3) * g(iy+3) * g(iz+3) &
                       + g(ix+4) * g(iy+4) * g(iz+4) &
                       + g(ix+5) * g(iy+5) * g(iz+5)
            end do
         case (7)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = g(ix) * g(iy) * g(iz) &
                       + g(ix+1) * g(iy+1) * g(iz+1) &
                       + g(ix+2) * g(iy+2) * g(iz+2) &
                       + g(ix+3) * g(iy+3) * g(iz+3) &
                       + g(ix+4) * g(iy+4) * g(iz+4) &
                       + g(ix+5) * g(iy+5) * g(iz+5) &
                       + g(ix+6) * g(iy+6) * g(iz+6)
            end do
         case (8)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = g(ix) * g(iy) * g(iz) &
                       + g(ix+1) * g(iy+1) * g(iz+1) &
                       + g(ix+2) * g(iy+2) * g(iz+2) &
                       + g(ix+3) * g(iy+3) * g(iz+3) &
                       + g(ix+4) * g(iy+4) * g(iz+4) &
                       + g(ix+5) * g(iy+5) * g(iz+5) &
                       + g(ix+6) * g(iy+6) * g(iz+6) &
                       + g(ix+7) * g(iy+7) * g(iz+7)
            end do
         case default
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               s = 0.0_dp
               do i = 0, nroots - 1
                  s = s + g(ix+i) * g(iy+i) * g(iz+i)
               end do
               gout(n) = s
            end do
         end select
      else
         select case (nroots)
         case (1)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = gout(n) + g(ix) * g(iy) * g(iz)
            end do
         case (2)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = gout(n) + (g(ix)   * g(iy)   * g(iz) &
                                  + g(ix+1) * g(iy+1) * g(iz+1))
            end do
         case (3)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = gout(n) + (g(ix)   * g(iy)   * g(iz) &
                                  + g(ix+1) * g(iy+1) * g(iz+1) &
                                  + g(ix+2) * g(iy+2) * g(iz+2))
            end do
         case (4)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = gout(n) + (g(ix)   * g(iy)   * g(iz) &
                                  + g(ix+1) * g(iy+1) * g(iz+1) &
                                  + g(ix+2) * g(iy+2) * g(iz+2) &
                                  + g(ix+3) * g(iy+3) * g(iz+3))
            end do
         case (5)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = gout(n) + (g(ix) * g(iy) * g(iz) &
                                  + g(ix+1) * g(iy+1) * g(iz+1) &
                                  + g(ix+2) * g(iy+2) * g(iz+2) &
                                  + g(ix+3) * g(iy+3) * g(iz+3) &
                                  + g(ix+4) * g(iy+4) * g(iz+4))
            end do
         case (6)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = gout(n) + (g(ix) * g(iy) * g(iz) &
                                  + g(ix+1) * g(iy+1) * g(iz+1) &
                                  + g(ix+2) * g(iy+2) * g(iz+2) &
                                  + g(ix+3) * g(iy+3) * g(iz+3) &
                                  + g(ix+4) * g(iy+4) * g(iz+4) &
                                  + g(ix+5) * g(iy+5) * g(iz+5))
            end do
         case (7)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = gout(n) + (g(ix) * g(iy) * g(iz) &
                                  + g(ix+1) * g(iy+1) * g(iz+1) &
                                  + g(ix+2) * g(iy+2) * g(iz+2) &
                                  + g(ix+3) * g(iy+3) * g(iz+3) &
                                  + g(ix+4) * g(iy+4) * g(iz+4) &
                                  + g(ix+5) * g(iy+5) * g(iz+5) &
                                  + g(ix+6) * g(iy+6) * g(iz+6))
            end do
         case (8)
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               gout(n) = gout(n) + (g(ix) * g(iy) * g(iz) &
                                  + g(ix+1) * g(iy+1) * g(iz+1) &
                                  + g(ix+2) * g(iy+2) * g(iz+2) &
                                  + g(ix+3) * g(iy+3) * g(iz+3) &
                                  + g(ix+4) * g(iy+4) * g(iz+4) &
                                  + g(ix+5) * g(iy+5) * g(iz+5) &
                                  + g(ix+6) * g(iy+6) * g(iz+6) &
                                  + g(ix+7) * g(iy+7) * g(iz+7))
            end do
         case default
            do n = 0, nf - 1
               ix = idx(n*3+0); iy = idx(n*3+1); iz = idx(n*3+2)
               s = 0.0_dp
               do i = 0, nroots - 1
                  s = s + g(ix+i) * g(iy+i) * g(iz+i)
               end do
               gout(n) = gout(n) + s
            end do
         end select
      end if
   end subroutine cint_gout2e

   function cint_2e_loop_nopt(gctr, gctroff, envs, ws, empty) result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer,  intent(in)    :: gctroff
      type(cint_env_vars), intent(inout) :: envs
      ! TARGET so the three loop-invariant blocks below can be bound once.
      ! Safe because nothing inside these loops can move ws%d: ws_ensure is
      ! called before any of it (see cint_workspace.f90) and ws_alloc_* only
      ! ever hands out space that is already there.
      type(cint_ws), intent(inout), target :: ws
      logical,  intent(inout) :: empty
      logical :: has_value

      integer :: i_sh, j_sh, k_sh, l_sh
      integer :: i_ctr, j_ctr, k_ctr, l_ctr
      integer :: i_prim, j_prim, k_prim, l_prim
      integer :: ai_o, aj_o, ak_o, al_o, ci_o, cj_o, ck_o, cl_o
      integer :: n_comp, nc, nf, lkl, lij
      integer :: leng, lenl, lenk, lenj, leni, len0
      integer :: og, ogout, ogctri, ogctrj, ogctrk, ogctrl
      integer :: oidx, olm, onon, osrt
      integer :: ip, jp, kp, lp, kij, kijb
      real(dp) :: expcutoff, rr_ij, rr_kl, akl, ekl, ccekl, log_rr_kl
      real(dp) :: eijcutoff, cutoff, expijkl, fac1i, fac1j, fac1k, fac1l
      real(dp) :: omega, omega2, theta, aij, dist_ij, rkl(0:2)
      ! The aliased flags.  The C keeps five in _empty[] and then adds a
      ! sixth by #define: gctrm is gctr, mempty is the CALLER'S empty, and
      ! m_ctr is n_comp.  So when n_comp == 1 the whole alias chain terminates
      ! on the caller's flag rather than on anything local -- which is how the
      ! single-component case reports that it produced something.  Missing
      ! that makes every contracted integral come back as if screened out.
      ! The optimizer, when the caller built one that fits.  It caches the
      ! logged coefficient bounds, the non-zero contraction map and the
      ! Cartesian index map -- three of the four things this prologue
      ! computes for every shell quartet.  All three are the same values by
      ! construction, so using them cannot move a result; the pair data is
      ! deliberately still computed here, because the C's cached copy is
      ! built with a different angular bound and would change which
      ! primitives survive screening.
      logical :: use_opt, use_f12
      ! L SHELLS.  `sp_x` is "this index carries s and p together"; `ctr_x`
      ! is "this index needs an explicit primitive-to-contracted step", which
      ! is the question the loop below actually asks and which an L shell
      ! answers yes to even at one contraction -- its coefficient is not a
      ! scalar and so cannot be folded into `fac`.  Every `x_ctr == 1` test
      ! the loop used to make is now `.not. ctr_x`, and with no L shell in
      ! the quartet the two are the same predicate on the same value.
      !
      ! `blk_x` is the stride at which that index's component changes within
      ! the (i,k,l,j) block -- see cint_prim_to_ctr_sp_0.
      logical :: sp_i, sp_j, sp_k, sp_l, any_sp
      logical :: ctr_i, ctr_j, ctr_k, ctr_l
      integer :: blk_i, blk_j, blk_k, blk_l
      integer :: okey, g0ok
      integer, parameter :: E_I = 0, E_J = 1, E_K = 2, E_L = 3, E_G = 4, E_M = 5
      logical :: flag(0:5)
      integer :: ii, ij, ik, il, ig, im
      ! The three workspace blocks the innermost loop passes down, bound once.
      ! `cint_g0_2e(ws%d(og:), ...)` builds an array descriptor at every call,
      ! and the innermost loop makes two and a half million of them; VTune put
      ! 0.03s of a 0.10s function on that one call statement.  The offsets are
      ! loop-invariant, so the descriptors can be too.
      ! CONTIGUOUS: true by construction -- all three are bound to whole
      ! slices of the workspace -- and it lets gfortran build the callee's
      ! descriptor without a stride check.
      real(dp), pointer, contiguous :: gp(:), goutp(:)
      integer,  pointer, contiguous :: idxp(:)
      ! THE EXPONENT AND COEFFICIENT ROWS, BOUND ONCE.  The C hoists these
      ! out of the loop nest -- `double *ai = env + bas(PTR_EXP,i_sh)` -- and
      ! then ai[ip] is one load.  Reaching through envs%env every time makes
      ! it three: the pointer component's descriptor, the offset, the value.
      ! Measured on the innermost line, envs%ai = envs%env(ai_o+ip) cost 30.4
      ! million instructions against the C's 5.0 million for the same
      ! assignment.
      real(dp), pointer, contiguous :: aip(:), ajp(:), akp(:), alp(:)
      real(dp), pointer, contiguous :: cip(:), cjp(:), ckp(:), clp(:)
      real(dp), pointer, contiguous :: lmkp(:), lmlp(:)
      type(pair_data), pointer, contiguous :: pdp(:)

      has_value = .false.
      i_sh = envs%shls(0); j_sh = envs%shls(1)
      k_sh = envs%shls(2); l_sh = envs%shls(3)
      i_ctr = envs%x_ctr(0); j_ctr = envs%x_ctr(1)
      k_ctr = envs%x_ctr(2); l_ctr = envs%x_ctr(3)
      i_prim = bas_of(envs, NPRIM_OF, i_sh)
      j_prim = bas_of(envs, NPRIM_OF, j_sh)
      k_prim = bas_of(envs, NPRIM_OF, k_sh)
      l_prim = bas_of(envs, NPRIM_OF, l_sh)
      ai_o = bas_of(envs, PTR_EXP, i_sh); aj_o = bas_of(envs, PTR_EXP, j_sh)
      ak_o = bas_of(envs, PTR_EXP, k_sh); al_o = bas_of(envs, PTR_EXP, l_sh)
      ci_o = bas_of(envs, PTR_COEFF, i_sh); cj_o = bas_of(envs, PTR_COEFF, j_sh)
      ck_o = bas_of(envs, PTR_COEFF, k_sh); cl_o = bas_of(envs, PTR_COEFF, l_sh)

      any_sp = envs%sp_mask /= 0
      sp_i = btest(envs%sp_mask, 0); sp_j = btest(envs%sp_mask, 1)
      sp_k = btest(envs%sp_mask, 2); sp_l = btest(envs%sp_mask, 3)
      ctr_i = i_ctr > 1 .or. sp_i; ctr_j = j_ctr > 1 .or. sp_j
      ctr_k = k_ctr > 1 .or. sp_k; ctr_l = l_ctr > 1 .or. sp_l
      blk_i = 1
      blk_k = envs%nfi
      blk_l = envs%nfi * envs%nfk
      blk_j = envs%nfi * envs%nfk * envs%nfl

      expcutoff = envs%expcutoff
      rr_ij = sum(envs%rirj * envs%rirj)
      rr_kl = sum(envs%rkrl * envs%rkrl)
      nf = envs%nf
      n_comp = envs%ncomp_e1 * envs%ncomp_e2 * envs%ncomp_tensor

      use_opt = .false.
      if (associated(ws%opt)) use_opt = cint_opt_usable(ws%opt, 4, envs)
      use_f12 = envs%g0_kind /= G0_COULOMB

      call ws_alloc_d(ws, i_prim + j_prim + k_prim + l_prim, olm)
      if (use_opt) then
         call ws_opt_log_maxc(ws, olm,          i_sh, i_prim)
         call ws_opt_log_maxc(ws, olm + i_prim, j_sh, j_prim)
      else
      ! merge(2,1,sp_*): an L shell's coefficient block is 2*nctr columns
      ! wide, and the bound this builds has to see the p ones too or a pair
      ! whose p contraction dominates gets screened on its s coefficient.
      call cint_log_max_pgto_coeff(ws%d(olm:), envs%env, ci_o, i_prim, &
                                   i_ctr*merge(2, 1, sp_i))
      call cint_log_max_pgto_coeff(ws%d(olm+i_prim:), envs%env, cj_o, j_prim, &
                                   j_ctr*merge(2, 1, sp_j))
      end if

      call ws_ensure_pd(ws, i_prim*j_prim)
      if (cint_set_pairdata(ws%p, envs%env, ai_o, envs%env, aj_o, envs%ri, envs%rj, &
                            ws%d(olm:), ws%d(olm+i_prim:), envs%li_ceil, envs%lj_ceil, &
                            i_prim, j_prim, rr_ij, expcutoff, &
                            envs%env(PTR_RANGE_OMEGA))) then
         return
      end if
      if (use_opt) then
         call ws_opt_log_maxc(ws, olm + i_prim + j_prim,          k_sh, k_prim)
         call ws_opt_log_maxc(ws, olm + i_prim + j_prim + k_prim, l_sh, l_prim)
      else
      call cint_log_max_pgto_coeff(ws%d(olm+i_prim+j_prim:), envs%env, ck_o, k_prim, &
                                   k_ctr*merge(2, 1, sp_k))
      call cint_log_max_pgto_coeff(ws%d(olm+i_prim+j_prim+k_prim:), envs%env, cl_o, l_prim, &
                                   l_ctr*merge(2, 1, sp_l))
      end if

      lkl = envs%lk_ceil + envs%ll_ceil
      akl = envs%env(ak_o + k_prim - 1) + envs%env(al_o + l_prim - 1)
      log_rr_kl = 1.7_dp - 1.5_dp * log(akl)
      omega = envs%env(PTR_RANGE_OMEGA)
      if (omega < 0.0_dp) then
         ! The screening estimate has to account for the attenuation, or the
         ! short-range integrals get cut too early.
         if (envs%rys_order > 1) then
            omega2 = omega * omega
            lij = envs%li_ceil + envs%lj_ceil
            if (lij > 0) then
               aij = envs%env(ai_o + i_prim - 1) + envs%env(aj_o + j_prim - 1)
               dist_ij = sqrt(rr_ij)
               theta = omega2 / (omega2 + aij)
               expcutoff = expcutoff + lij * &
                  log((dist_ij + theta*R_GUESS + 1.0_dp) / (dist_ij + 1.0_dp))
            end if
            if (lkl > 0) then
               theta = omega2 / (omega2 + akl)
               log_rr_kl = log_rr_kl + lkl * log(sqrt(rr_kl) + theta*R_GUESS + 1.0_dp)
            end if
         end if
      else
         if (lkl > 0) log_rr_kl = log_rr_kl + lkl * log(sqrt(rr_kl) + 1.0_dp)
      end if

      okey = -1
      if (use_opt .and. .not. any_sp) then
         ! NOT FOR L SHELLS.  The cache is keyed on the four angular momenta
         ! and built from a fake basis that has one shell per l, so it has no
         ! way to hold a second map for the same l with an s prepended.  The
         ! map is built here instead, which is what the unoptimised path does
         ! for every quartet anyway.
         okey = ws%opt%idx_off(cint_opt_idx_key(envs, 4))
         ! l above the l_allow gen_idx was built with leaves a hole in the
         ! table.  The C indexes into it anyway and dereferences NULL; here
         ! the hole says OPT_NOVALUE and the map is built as usual.
      end if
      if (okey == OPT_NOVALUE) okey = -1
      if (okey < 0) then
         call ws_alloc_i(ws, nf*3, oidx)
         call cint_g2e_index_xyz(ws%i(oidx:), envs)
      end if

      call ws_alloc_i(ws, i_prim + j_prim + k_prim + l_prim, onon)
      call ws_alloc_i(ws, i_prim*i_ctr + j_prim*j_ctr + k_prim*k_ctr + l_prim*l_ctr, osrt)
      if (use_opt) then
         call ws_opt_non0(ws, onon, osrt, i_sh, i_prim, i_ctr, 0, 0)
         call ws_opt_non0(ws, onon, osrt, j_sh, j_prim, j_ctr, i_prim, i_prim*i_ctr)
         call ws_opt_non0(ws, onon, osrt, k_sh, k_prim, k_ctr, i_prim+j_prim, &
                            i_prim*i_ctr+j_prim*j_ctr)
         call ws_opt_non0(ws, onon, osrt, l_sh, l_prim, l_ctr, i_prim+j_prim+k_prim, &
                            i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr)
      else
      call cint_non0coeff_byshell(ws%i(osrt:), ws%i(onon:), envs%env, ci_o, i_prim, i_ctr)
      call cint_non0coeff_byshell(ws%i(osrt+i_prim*i_ctr:), ws%i(onon+i_prim:), &
                                  envs%env, cj_o, j_prim, j_ctr)
      call cint_non0coeff_byshell(ws%i(osrt+i_prim*i_ctr+j_prim*j_ctr:), &
                                  ws%i(onon+i_prim+j_prim:), envs%env, ck_o, k_prim, k_ctr)
      call cint_non0coeff_byshell(ws%i(osrt+i_prim*i_ctr+j_prim*j_ctr+k_prim*k_ctr:), &
                                  ws%i(onon+i_prim+j_prim+k_prim:), envs%env, cl_o, l_prim, l_ctr)
      end if

      nc = i_ctr * j_ctr * k_ctr * l_ctr
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenl = nf * nc * n_comp
      lenk = nf * i_ctr * j_ctr * k_ctr * n_comp
      lenj = nf * i_ctr * j_ctr * n_comp
      leni = nf * i_ctr * n_comp
      len0 = nf * n_comp

      call ws_alloc_d(ws, leng, og)
      ! ALIAS_ADDR_IF_EQUAL: each level shares its successor's buffer when the
      ! contraction count is 1, and shares its empty flag with it too.
      flag = .true.
      flag(E_M) = empty
      ii = E_I; ij = E_J; ik = E_K; il = E_L; ig = E_G; im = E_M
      if (n_comp == 1) then
         ogctrl = gctroff; il = im
      else
         call ws_alloc_d(ws, lenl, ogctrl)
      end if
      ! ctr_x, not x_ctr == 1: a level that has a contraction step to do
      ! cannot share its successor's buffer, and an L shell has one whatever
      ! its contraction count.
      if (.not. ctr_l) then
         ogctrk = ogctrl; ik = il
      else
         call ws_alloc_d(ws, lenk, ogctrk)
      end if
      if (.not. ctr_k) then
         ogctrj = ogctrk; ij = ik
      else
         call ws_alloc_d(ws, lenj, ogctrj)
      end if
      if (.not. ctr_j) then
         ogctri = ogctrj; ii = ij
      else
         call ws_alloc_d(ws, leni, ogctri)
      end if
      if (.not. ctr_i) then
         ogout = ogctri; ig = ii
      else
         call ws_alloc_d(ws, len0, ogout)
      end if

      aip(0:) => envs%env(ai_o:); ajp(0:) => envs%env(aj_o:)
      akp(0:) => envs%env(ak_o:); alp(0:) => envs%env(al_o:)
      cip(0:) => envs%env(ci_o:); cjp(0:) => envs%env(cj_o:)
      ckp(0:) => envs%env(ck_o:); clp(0:) => envs%env(cl_o:)
      lmkp(0:) => ws%d(olm+i_prim+j_prim:)
      lmlp(0:) => ws%d(olm+i_prim+j_prim+k_prim:)
      pdp(0:) => ws%p

      gp(0:)    => ws%d(og:og+leng-1)
      goutp(0:) => ws%d(ogout:ogout+len0-1)
      if (okey >= 0) then
         idxp(0:) => ws%opt%idx(okey:okey+nf*3-1)
      else
         idxp(0:) => ws%i(oidx:oidx+nf*3-1)
      end if

      if (.not. (ctr_i .or. ctr_j .or. ctr_k .or. ctr_l)) then
      ! Segmented shells: the C's CINT2e_1111_loop.
      !
      ! When every contraction count is one, the whole apparatus below
      ! collapses.  The five `ogctr*` buffers alias onto one another and onto
      ! `gctr` (see the chain above, each `== 1` branch folding the next), so
      ! `f_gout` already writes where the answer belongs; the six `flag`
      ! entries alias onto one -- `ig`, whichever that turned out to be -- so
      ! there is one empty flag rather than six; and every `prim_to_ctr` is
      ! guarded by a `> 1` that is false.
      !
      ! The test is `ctr_*` rather than `x_ctr == 1` for the L shells: one of
      ! those has a real contraction step at every level whatever its
      ! contraction count, because its coefficient is a vector over the four
      ! components and cannot be folded into `fac`.  With no L shell in the
      ! quartet the two predicates are the same.
      !
      ! `n_comp` is deliberately not part of the test.  It only decides which
      ! flag the chain lands on and whether the transpose after the loop runs;
      ! the aliasing is driven by the four contraction counts alone.  Writing
      ! `flag(ig)` rather than naming a flag keeps that true, and the
      ! transpose below is left to the guard it already has -- which is what
      ! lets this cover the derivative integrals, where n_comp is 3 and the
      ! gradient spends most of its time.
      !
      ! The generic nest still executes all of it -- four branches per
      ! primitive level, an aliased-index store per surviving quartet, and
      ! four epilogue tests per level -- to reach that conclusion every time.
      ! On a Pople basis it reaches it for every shell in the molecule, since
      ! 6-31G is segmented throughout.
      !
      ! Bit-identical by construction rather than by measurement: this is the
      ! same products in the same order, with the branches that cannot be
      ! taken removed.  `catalogue_check` still holds it to the last bit.
      !
      ! Two loops, not one, is what libcint does here too, and for the same
      ! reason: the general one has to stay for generally-contracted bases,
      ! where the contraction is the whole point.  Deleting it is what makes
      ! qcint slower than libcint on exactly this case.
      do lp = 0, l_prim - 1
         envs%al = alp(lp)
         fac1l = envs%common_factor * clp(lp)

         do kp = 0, k_prim - 1
            akl = akp(kp) + alp(lp)
            ekl = rr_kl * akp(kp) * alp(lp) / akl
            ccekl = ekl - log_rr_kl - lmkp(kp) - lmlp(lp)
            if (ccekl > expcutoff) cycle
            envs%ak = akp(kp)
            rkl(0) = (akp(kp)*envs%rk(0) + alp(lp)*envs%rl(0)) / akl
            rkl(1) = (akp(kp)*envs%rk(1) + alp(lp)*envs%rl(1)) / akl
            rkl(2) = (akp(kp)*envs%rk(2) + alp(lp)*envs%rl(2)) / akl
            eijcutoff = expcutoff - ccekl
            ekl = exp(-ekl)
            fac1k = fac1l * ckp(kp)

            do jp = 0, j_prim - 1
               envs%aj = ajp(jp)
               fac1j = fac1k * cjp(jp)

               kijb = jp * i_prim
               do ip = 0, i_prim - 1
                  kij = kijb + ip
                  if (pdp(kij)%cceij > eijcutoff) cycle
                  envs%ai = aip(ip)
                  cutoff = eijcutoff - pdp(kij)%cceij
                  expijkl = pdp(kij)%eij * ekl
                  fac1i = fac1j * cip(ip) * expijkl
                  envs%fac = fac1i
                  if (use_f12) then
                     g0ok = f12_hook(gp, pdp(kij)%rij, rkl, cutoff, envs)
                  else
                     g0ok = cint_g0_2e(gp, pdp(kij)%rij, rkl, cutoff, envs)
                  end if
                  if (g0ok /= 0) then
                     call envs%f_gout(goutp, gp, idxp, envs, merge(1, 0, flag(ig)))
                     flag(ig) = .false.
                  end if
               end do
            end do
         end do
      end do
      else
      do lp = 0, l_prim - 1
         envs%al = alp(lp)
         if (.not. ctr_l) then
            fac1l = envs%common_factor * clp(lp)
         else
            fac1l = envs%common_factor
            flag(ik) = .true.
         end if

         do kp = 0, k_prim - 1
            akl = akp(kp) + alp(lp)
            ekl = rr_kl * akp(kp) * alp(lp) / akl
            ccekl = ekl - log_rr_kl - lmkp(kp) - lmlp(lp)
            if (ccekl > expcutoff) cycle
            envs%ak = akp(kp)
            rkl(0) = (akp(kp)*envs%rk(0) + alp(lp)*envs%rl(0)) / akl
            rkl(1) = (akp(kp)*envs%rk(1) + alp(lp)*envs%rl(1)) / akl
            rkl(2) = (akp(kp)*envs%rk(2) + alp(lp)*envs%rl(2)) / akl
            eijcutoff = expcutoff - ccekl
            ekl = exp(-ekl)
            if (.not. ctr_k) then
               fac1k = fac1l * ckp(kp)
            else
               fac1k = fac1l
               flag(ij) = .true.
            end if

            do jp = 0, j_prim - 1
               envs%aj = ajp(jp)
               if (.not. ctr_j) then
                  fac1j = fac1k * cjp(jp)
               else
                  fac1j = fac1k
                  flag(ii) = .true.
               end if

               kijb = jp * i_prim
               do ip = 0, i_prim - 1
                  ! One add, not an increment on each of two paths.  The C
                  ! bumps pdata_ij in the loop header, so the screened path
                  ! and the surviving one share it.
                  kij = kijb + ip
                  if (pdp(kij)%cceij > eijcutoff) cycle
                  envs%ai = aip(ip)
                  cutoff = eijcutoff - pdp(kij)%cceij
                  expijkl = pdp(kij)%eij * ekl
                  if (.not. ctr_i) then
                     fac1i = fac1j * cip(ip) * expijkl
                  else
                     fac1i = fac1j * expijkl
                  end if
                  envs%fac = fac1i
                  ! A branch that goes the same way for a whole run, against
                  ! an indirect call the C takes unconditionally.
                  if (use_f12) then
                     g0ok = f12_hook(gp, pdp(kij)%rij, rkl, cutoff, envs)
                  else
                     g0ok = cint_g0_2e(gp, pdp(kij)%rij, rkl, cutoff, envs)
                  end if
                  if (g0ok /= 0) then
                     call envs%f_gout(goutp, gp, idxp, envs, merge(1, 0, flag(ig)))
                     if (ctr_i) then
                        if (sp_i) then
                           if (flag(ii)) then
                              call cint_prim_to_ctr_sp_0(ws%d, ogctri, ogout, envs%env, ci_o+ip, &
                                                         len0, blk_i, i_prim, i_ctr)
                           else
                              call cint_prim_to_ctr_sp_1(ws%d, ogctri, ogout, envs%env, ci_o+ip, &
                                                         len0, blk_i, i_prim, i_ctr)
                           end if
                        else if (flag(ii)) then
                           call cint_prim_to_ctr_0(ws%d, ogctri, ogout, envs%env, ci_o+ip, &
                                                   len0, i_prim, i_ctr)
                        else
                           call cint_prim_to_ctr_1(ws%d, ogctri, ogout, envs%env, ci_o+ip, &
                                                   len0, i_prim, i_ctr, ws%i(onon+ip), ws%i, osrt + ip*i_ctr)
                        end if
                     end if
                     flag(ii) = .false.
                  end if
               end do

               if (.not. flag(ii)) then
                  if (ctr_j) then
                     if (sp_j) then
                        if (flag(ij)) then
                           call cint_prim_to_ctr_sp_0(ws%d, ogctrj, ogctri, envs%env, cj_o+jp, &
                                                      leni, blk_j, j_prim, j_ctr)
                        else
                           call cint_prim_to_ctr_sp_1(ws%d, ogctrj, ogctri, envs%env, cj_o+jp, &
                                                      leni, blk_j, j_prim, j_ctr)
                        end if
                     else if (flag(ij)) then
                        call cint_prim_to_ctr_0(ws%d, ogctrj, ogctri, envs%env, cj_o+jp, &
                                                leni, j_prim, j_ctr)
                     else
                        call cint_prim_to_ctr_1(ws%d, ogctrj, ogctri, envs%env, cj_o+jp, &
                                                leni, j_prim, j_ctr, ws%i(onon+i_prim+jp), ws%i, osrt + i_prim*i_ctr + jp*j_ctr)
                     end if
                  end if
                  flag(ij) = .false.
               end if
            end do

            if (.not. flag(ij)) then
               if (ctr_k) then
                  if (sp_k) then
                     if (flag(ik)) then
                        call cint_prim_to_ctr_sp_0(ws%d, ogctrk, ogctrj, envs%env, ck_o+kp, &
                                                   lenj, blk_k, k_prim, k_ctr)
                     else
                        call cint_prim_to_ctr_sp_1(ws%d, ogctrk, ogctrj, envs%env, ck_o+kp, &
                                                   lenj, blk_k, k_prim, k_ctr)
                     end if
                  else if (flag(ik)) then
                     call cint_prim_to_ctr_0(ws%d, ogctrk, ogctrj, envs%env, ck_o+kp, &
                                             lenj, k_prim, k_ctr)
                  else
                     call cint_prim_to_ctr_1(ws%d, ogctrk, ogctrj, envs%env, ck_o+kp, &
                                             lenj, k_prim, k_ctr, ws%i(onon+i_prim+j_prim+kp), ws%i, &
                                                                       osrt + i_prim*i_ctr + j_prim*j_ctr + kp*k_ctr)
                  end if
               end if
               flag(ik) = .false.
            end if
         end do

         if (.not. flag(ik)) then
            if (ctr_l) then
               if (sp_l) then
                  if (flag(il)) then
                     call cint_prim_to_ctr_sp_0(ws%d, ogctrl, ogctrk, envs%env, cl_o+lp, &
                                                lenk, blk_l, l_prim, l_ctr)
                  else
                     call cint_prim_to_ctr_sp_1(ws%d, ogctrl, ogctrk, envs%env, cl_o+lp, &
                                                lenk, blk_l, l_prim, l_ctr)
                  end if
               else if (flag(il)) then
                  call cint_prim_to_ctr_0(ws%d, ogctrl, ogctrk, envs%env, cl_o+lp, &
                                          lenk, l_prim, l_ctr)
               else
                  call cint_prim_to_ctr_1(ws%d, ogctrl, ogctrk, envs%env, cl_o+lp, &
                                          lenk, l_prim, l_ctr, ws%i(onon+i_prim+j_prim+k_prim+lp), ws%i, &
                                                                    osrt + i_prim*i_ctr + j_prim*j_ctr + k_prim*k_ctr + lp*l_ctr)
               end if
            end if
            flag(il) = .false.
         end if
      end do
      end if

      if (n_comp > 1 .and. .not. flag(il)) then
         if (flag(im)) then
            call cint_dmat_transpose(gctr(gctroff:), ws%d(ogctrl:), nf*nc, n_comp)
            flag(im) = .false.
         else
            call cint_dplus_transpose(gctr(gctroff:), ws%d(ogctrl:), nf*nc, n_comp)
         end if
      end if
      empty = flag(im)
      has_value = .not. empty
   end function cint_2e_loop_nopt

   ! PRIM2CTR is a macro in the C (src/cint2e.c), and the four call sites
   ! below expand it the same way.  It was a subroutine here, and the wrapper
   ! alone -- prologue, the empty-flag test, epilogue -- cost 32 million
   ! instructions across 350,616 calls, none of which the C pays.

   ! How many spherical functions one index of the quartet contributes.
   !
   ! An L shell contributes four -- one s and three p -- and they ARE the
   ! Cartesian four, in the same order, because cint_c2s_bra_sph and
   ! cint_c2s_ket_sph both return RESULT_IN_GCART below l = 2.  So an L shell
   ! needs no cart-to-spherical transform at all, and the only thing the four
   ! transform stages below have to be told about one is its dimension.
   pure integer function sph_dim(envs, x, l) result(d)
      type(cint_env_vars), intent(in) :: envs
      integer,             intent(in) :: x, l
      if (btest(envs%sp_mask, x)) then
         d = NF_SP
      else
         d = l*2 + 1
      end if
   end function sph_dim

   pure subroutine int2e_cache_size(envs, nd, ni)
      type(cint_env_vars), intent(in) :: envs
      integer, intent(out) :: nd, ni
      integer :: ip, jp, kp, lp, ic, jc, kc, lc, nc, n_comp, nf
      integer :: leng, lenl, lenk, lenj, leni, len0, buflen

      ip = bas_of(envs, NPRIM_OF, envs%shls(0))
      jp = bas_of(envs, NPRIM_OF, envs%shls(1))
      kp = bas_of(envs, NPRIM_OF, envs%shls(2))
      lp = bas_of(envs, NPRIM_OF, envs%shls(3))
      ic = envs%x_ctr(0); jc = envs%x_ctr(1)
      kc = envs%x_ctr(2); lc = envs%x_ctr(3)
      nf = envs%nf
      nc = ic*jc*kc*lc
      n_comp = envs%ncomp_e1 * envs%ncomp_e2 * envs%ncomp_tensor

      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenl = nf * nc * n_comp
      lenk = nf * ic*jc*kc * n_comp
      lenj = nf * ic*jc * n_comp
      leni = nf * ic * n_comp
      len0 = nf * n_comp
      ! the four c2s_sph_2e1 scratch blocks
      buflen = envs%nfi * envs%nfk * envs%nfl * sph_dim(envs, 1, envs%j_l)

      nd = nf*nc*n_comp + ip+jp+kp+lp &
         + leng + lenl + lenk + lenj + leni + len0 &
         + 4*buflen + 64
      ni = nf*3 + ip+jp+kp+lp + ip*ic + jp*jc + kp*kc + lp*lc + 64
   end subroutine int2e_cache_size

   function cint_2e_drv(out, dims, envs, ws, c2s_kind) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      integer,  intent(in)    :: c2s_kind
      logical :: has_value
      integer :: dmark, imark

      integer :: nc, n_comp, ogctr, n, nout, counts(0:3), nd, ni
      logical :: empty, dummy

      nc = envs%nf * envs%x_ctr(0) * envs%x_ctr(1) * envs%x_ctr(2) * envs%x_ctr(3)
      n_comp = envs%ncomp_e1 * envs%ncomp_e2 * envs%ncomp_tensor

      call int2e_cache_size(envs, nd, ni)
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*n_comp, ogctr)
      ! Not zeroed.  CINT2e_drv does not zero its gctr either: the loop's
      ! empty-flag protocol makes the first write to every component an
      ! assignment rather than an accumulation, and anything never written
      ! is never read, because has_value says so.  Zeroing it was 37.3
      ! million instructions -- one memset per shell quartet, 1% of the
      ! port's total -- spent restating what the flags already guarantee.

      empty = .true.
      dummy = cint_2e_loop_nopt(ws%d, ogctr, envs, ws, empty)
      has_value = .not. empty

      select case (c2s_kind)
      case (C2S_SPH_2E1)
         counts(0) = sph_dim(envs, 0, envs%i_l) * envs%x_ctr(0)
         counts(1) = sph_dim(envs, 1, envs%j_l) * envs%x_ctr(1)
         counts(2) = sph_dim(envs, 2, envs%k_l) * envs%x_ctr(2)
         counts(3) = sph_dim(envs, 3, envs%l_l) * envs%x_ctr(3)
      case default
         counts(0) = envs%nfi * envs%x_ctr(0)
         counts(1) = envs%nfj * envs%x_ctr(1)
         counts(2) = envs%nfk * envs%x_ctr(2)
         counts(3) = envs%nfl * envs%x_ctr(3)
      end select
      nout = dims(0) * dims(1) * dims(2) * dims(3)

      if (has_value) then
         call ws_mark(ws, dmark, imark)
         do n = 0, n_comp - 1
            call ws_rewind(ws, dmark, imark)
            select case (c2s_kind)
            case (C2S_CART_2E1)
               call apply_c2s_cart_2e1(out(nout*n:), ws%d(ogctr+nc*n:), dims, envs)
            case (C2S_SPH_2E1)
               call apply_c2s_sph_2e1(out(nout*n:), ogctr+nc*n, dims, envs, ws)
            case default
               error stop "cint_2e_drv: spinor c2s is D9"
            end select
         end do
      else
         do n = 0, n_comp - 1
            call apply_c2s_dset0_2e(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_2e_drv

   ! Reorder the (i,k,l,j) block the g array produces into the caller's
   ! (i,j,k,l) layout.
   pure subroutine dcopy_iklj(fijkl, foff, gctr, goff, ni, nj, nk, nl, mi, mj, mk, ml)
      real(dp), intent(inout) :: fijkl(0:)
      real(dp), intent(in)    :: gctr(0:)
      integer,  intent(in)    :: foff, goff, ni, nj, nk, nl, mi, mj, mk, ml
      integer :: i, j, k, l, nij, nijk, mik, mikl, fb, gb
      nij = ni*nj; nijk = nij*nk; mik = mi*mk; mikl = mik*ml
      fb = foff; gb = goff
      ! UNROLLED ON mi, as the C is (src/cart2sph.c:4570).  mi is the bra
      ! dimension: 1, 3, 5, 7 for s, p, d, f spherical and 6 for Cartesian d,
      ! so the general loop is a two-trip loop almost every time -- and
      ! gfortran turned it into a memcpy call, 1,656,493 of them for 20.2
      ! million instructions of call overhead moving a handful of doubles.
      select case (mi)
      case (1)
         do l = 0, ml - 1
            do k = 0, mk - 1
               do j = 0, mj - 1
                  fijkl(fb + k*nij + ni*j) = gctr(gb + k*mi + mikl*j)
               end do
            end do
            fb = fb + nijk
            gb = gb + mik
         end do
      case (3)
         do l = 0, ml - 1
            do k = 0, mk - 1
               do j = 0, mj - 1
                  fijkl(fb + k*nij + ni*j + 0) = gctr(gb + k*mi + mikl*j + 0)
                  fijkl(fb + k*nij + ni*j + 1) = gctr(gb + k*mi + mikl*j + 1)
                  fijkl(fb + k*nij + ni*j + 2) = gctr(gb + k*mi + mikl*j + 2)
               end do
            end do
            fb = fb + nijk
            gb = gb + mik
         end do
      case (4)
         ! An L shell's bra dimension, which the C has no case for because it
         ! cannot express one.  The general loop below produces it correctly
         ! and gfortran turns it into a memcpy call, for the same handful of
         ! doubles the other cases are unrolled to avoid.
         do l = 0, ml - 1
            do k = 0, mk - 1
               do j = 0, mj - 1
                  fijkl(fb + k*nij + ni*j + 0) = gctr(gb + k*mi + mikl*j + 0)
                  fijkl(fb + k*nij + ni*j + 1) = gctr(gb + k*mi + mikl*j + 1)
                  fijkl(fb + k*nij + ni*j + 2) = gctr(gb + k*mi + mikl*j + 2)
                  fijkl(fb + k*nij + ni*j + 3) = gctr(gb + k*mi + mikl*j + 3)
               end do
            end do
            fb = fb + nijk
            gb = gb + mik
         end do
      case (5)
         do l = 0, ml - 1
            do k = 0, mk - 1
               do j = 0, mj - 1
                  fijkl(fb + k*nij + ni*j + 0) = gctr(gb + k*mi + mikl*j + 0)
                  fijkl(fb + k*nij + ni*j + 1) = gctr(gb + k*mi + mikl*j + 1)
                  fijkl(fb + k*nij + ni*j + 2) = gctr(gb + k*mi + mikl*j + 2)
                  fijkl(fb + k*nij + ni*j + 3) = gctr(gb + k*mi + mikl*j + 3)
                  fijkl(fb + k*nij + ni*j + 4) = gctr(gb + k*mi + mikl*j + 4)
               end do
            end do
            fb = fb + nijk
            gb = gb + mik
         end do
      case (6)
         do l = 0, ml - 1
            do k = 0, mk - 1
               do j = 0, mj - 1
                  fijkl(fb + k*nij + ni*j + 0) = gctr(gb + k*mi + mikl*j + 0)
                  fijkl(fb + k*nij + ni*j + 1) = gctr(gb + k*mi + mikl*j + 1)
                  fijkl(fb + k*nij + ni*j + 2) = gctr(gb + k*mi + mikl*j + 2)
                  fijkl(fb + k*nij + ni*j + 3) = gctr(gb + k*mi + mikl*j + 3)
                  fijkl(fb + k*nij + ni*j + 4) = gctr(gb + k*mi + mikl*j + 4)
                  fijkl(fb + k*nij + ni*j + 5) = gctr(gb + k*mi + mikl*j + 5)
               end do
            end do
            fb = fb + nijk
            gb = gb + mik
         end do
      case (7)
         do l = 0, ml - 1
            do k = 0, mk - 1
               do j = 0, mj - 1
                  fijkl(fb + k*nij + ni*j + 0) = gctr(gb + k*mi + mikl*j + 0)
                  fijkl(fb + k*nij + ni*j + 1) = gctr(gb + k*mi + mikl*j + 1)
                  fijkl(fb + k*nij + ni*j + 2) = gctr(gb + k*mi + mikl*j + 2)
                  fijkl(fb + k*nij + ni*j + 3) = gctr(gb + k*mi + mikl*j + 3)
                  fijkl(fb + k*nij + ni*j + 4) = gctr(gb + k*mi + mikl*j + 4)
                  fijkl(fb + k*nij + ni*j + 5) = gctr(gb + k*mi + mikl*j + 5)
                  fijkl(fb + k*nij + ni*j + 6) = gctr(gb + k*mi + mikl*j + 6)
               end do
            end do
            fb = fb + nijk
            gb = gb + mik
         end do
      case default
         do l = 0, ml - 1
            do k = 0, mk - 1
               do j = 0, mj - 1
                  do i = 0, mi - 1
                     fijkl(fb + k*nij + ni*j + i) = gctr(gb + k*mi + mikl*j + i)
                  end do
               end do
            end do
            fb = fb + nijk
            gb = gb + mik
         end do
      end select
   end subroutine dcopy_iklj

   subroutine apply_c2s_cart_2e1(fijkl, gctr, dims, envs)
      real(dp), intent(inout) :: fijkl(0:)
      real(dp), intent(in)    :: gctr(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(in) :: envs
      integer :: ic, jc, kc, lc, ni, nj, nk, ofj, ofk, ofl, base, gb
      ni = dims(0); nj = dims(1); nk = dims(2)
      ofj = ni * envs%nfj
      ofk = ni * nj * envs%nfk
      ofl = ni * nj * nk * envs%nfl
      gb = 0
      do lc = 0, envs%x_ctr(3) - 1
      do kc = 0, envs%x_ctr(2) - 1
      do jc = 0, envs%x_ctr(1) - 1
      do ic = 0, envs%x_ctr(0) - 1
         base = ofl*lc + ofk*kc + ofj*jc + envs%nfi*ic
         call dcopy_iklj(fijkl, base, gctr, gb, ni, nj, nk, dims(3), &
                         envs%nfi, envs%nfj, envs%nfk, envs%nfl)
         gb = gb + envs%nf
      end do
      end do
      end do
      end do
   end subroutine apply_c2s_cart_2e1

   ! ONE ARRAY, TWO OFFSETS.  Both the source and the destination are the
   ! workspace, and passing the whole of it as an intent(in) dummy and again
   ! as an intent(inout) one is the aliasing violation this port keeps
   ! rediscovering.  The slices do not overlap; saying so with offsets is
   ! what makes that a fact rather than a hope.
   subroutine sph2e_inner(buf, soff, coff, l, nbra, ncall, sizsph, sizcart, loc)
      real(dp), intent(inout) :: buf(0:)
      integer,  intent(in)    :: soff, coff, l, nbra, ncall, sizsph, sizcart
      integer,  intent(out)   :: loc
      integer :: n, r
      if (l <= 1) then
         loc = RESULT_IN_GCART
         return
      end if
      do n = 0, ncall - 1
         r = cint_c2s_ket_sph(buf(soff + n*sizsph:), buf(coff + n*sizcart:), &
                              nbra, nbra, l)
      end do
      loc = RESULT_IN_GSPH
   end subroutine sph2e_inner

   ! gctr arrives as an offset into the workspace rather than a section of
   ! it, because sph2e_inner below indexes the same buffer and the two have
   ! to share a base.
   subroutine apply_c2s_sph_2e1(out, goff, dims, envs, ws)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: goff, dims(0:)
      type(cint_env_vars), intent(in) :: envs
      type(cint_ws), intent(inout) :: ws

      integer :: i_l, j_l, k_l, l_l, di, dj, dk, dl

      integer :: icp
      integer :: ni, nj, nk, nfi, nfk, nfl, nfik, nfikl, dlj
      integer :: ofj, ofk, ofl, ic, jc, kc, lc, buflen
      integer :: o1, o2, o3, o4, gb, base, cur, loc

      i_l = envs%i_l; j_l = envs%j_l; k_l = envs%k_l; l_l = envs%l_l
      di = sph_dim(envs, 0, i_l); dj = sph_dim(envs, 1, j_l)
      dk = sph_dim(envs, 2, k_l); dl = sph_dim(envs, 3, l_l)
      ni = dims(0); nj = dims(1); nk = dims(2)
      nfi = envs%nfi; nfk = envs%nfk; nfl = envs%nfl
      nfik = nfi*nfk; nfikl = nfik*nfl; dlj = dl*dj
      ofj = ni*dj; ofk = ni*nj*dk; ofl = ni*nj*nk*dl
      buflen = nfikl*dj
      call ws_alloc_d(ws, buflen, o1)
      call ws_alloc_d(ws, buflen, o2)
      call ws_alloc_d(ws, buflen, o3)
      call ws_alloc_d(ws, buflen, o4)

      gb = goff
      do lc = 0, envs%x_ctr(3) - 1
      do kc = 0, envs%x_ctr(2) - 1
      do jc = 0, envs%x_ctr(1) - 1
      do ic = 0, envs%x_ctr(0) - 1
         ! j, then l, then k, then i -- the order the C uses, which is the
         ! order the (i,k,l,j) layout makes cheapest.
         ! THREAD THE OFFSET, DO NOT NORMALISE THE BUFFER.  The C carries a
         ! tmp1 pointer through the four stages: each returns either its own
         ! output buffer or its input untouched, and the next stage reads
         ! whichever came back (src/cart2sph.c, c2s_sph_2e1).  Copying to put
         ! the result where the next stage expected it cost two loops here,
         ! and before that a malloc and a memcpy for each -- the section
         ! assignment they replaced had ws%d on both sides, so gfortran could
         ! not prove the ranges disjoint and allocated a temporary.
         loc = cint_c2s_ket_sph(ws%d(o1:), ws%d(gb:), nfikl, nfikl, j_l)
         if (loc == RESULT_IN_GCART) then
            cur = gb
         else
            cur = o1
         end if

         call sph2e_inner(ws%d, o2, cur, l_l, nfik, dj, nfik*dl, nfikl, loc)
         if (loc /= RESULT_IN_GCART) cur = o2

         call sph2e_inner(ws%d, o3, cur, k_l, nfi, dlj, nfi*dk, nfik, loc)
         if (loc /= RESULT_IN_GCART) cur = o3

         loc = cint_c2s_bra_sph(ws%d(o4:), dk*dlj, ws%d(cur:), i_l)
         if (loc == RESULT_IN_GCART) then
            base = ofl*lc + ofk*kc + ofj*jc + di*ic
            call dcopy_iklj(out, base, ws%d, cur, ni, nj, nk, dims(3), di, dj, dk, dl)
         else
            base = ofl*lc + ofk*kc + ofj*jc + di*ic
            call dcopy_iklj(out, base, ws%d, o4, ni, nj, nk, dims(3), di, dj, dk, dl)
         end if
         gb = gb + envs%nf
      end do
      end do
      end do
      end do
   end subroutine apply_c2s_sph_2e1

   pure subroutine apply_c2s_dset0_2e(out, dims, counts)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), counts(0:)
      integer :: i, j, k, l, ni, nij, nijk
      ni = dims(0); nij = ni*dims(1); nijk = nij*dims(2)
      do l = 0, counts(3) - 1
      do k = 0, counts(2) - 1
      do j = 0, counts(1) - 1
      do i = 0, counts(0) - 1
         out(l*nijk + k*nij + j*ni + i) = 0.0_dp
      end do
      end do
      end do
      end do
   end subroutine apply_c2s_dset0_2e

   ! ---- entry points -------------------------------------------------

   function int2e_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]

      call cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_CART_2E1)
   end function int2e_cart

   function int2e_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]

      call cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      has_value = cint_2e_drv(out, dims, envs, ws, C2S_SPH_2E1)
   end function int2e_sph

   ! ---- optimizers ------------------------------------------------------

   subroutine int2e_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_all_2e_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int2e_optimizer


end module cint_2e
