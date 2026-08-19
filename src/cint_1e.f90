!
! The one-electron driver, and int1e_ovlp on top of it.
!
! Ported from src/cint1e.c.  This is the layer that turns a shell pair into a
! block of integrals: loop the primitive pairs, build the G tensor for each,
! contract it into the output, and finally transform Cartesian to spherical.
!
! Three things differ from the C, all noted where they occur:
!
!   * scratch comes from a cint_ws workspace rather than MALLOC_INSTACK over a
!     caller-supplied cache (see cint_workspace.f90);
!   * the c2s routine is chosen by an enum rather than by comparing function
!     pointers, which is what CINT1e_drv does to work out its own output
!     dimensions;
!   * the output is real or complex depending on that enum.  Only the real
!     arms exist so far -- the spinor ones are D9 -- but the signature carries
!     both from the start, because retrofitting a complex return through this
!     layer later is exactly the rework PORT_TO_FORTRAN.md 4.5 warns about.
!
module cint_1e
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws, ws_ensure, ws_ensure_pd, ws_alloc_d, ws_alloc_i, ws_mark, ws_rewind, &
                             ws_opt_log_maxc, ws_opt_non0, ws_opt_idx
   use cint_screen,    only: pair_data, cint_set_pairdata, &
                             cint_log_max_pgto_coeff, cint_non0coeff_byshell
   use cint_g1e
   ! The optimizer builder comes in with the bare `use cint_g1e` above; naming
   ! it again here adds no visibility and ifx says so, remark #6536.
   use cint_bas,       only: cint_square_dist, &
                             NPRIM_OF, PTR_EXP, PTR_COEFF
   use cint_blas,      only: cint_dmat_transpose
   use cint_cart2sph,  only: cint_c2s_bra_sph, cint_c2s_ket_sph, RESULT_IN_GCART
   implicit none
   private

   public :: cint_gout1e, cint_gout1e_nuc, cint_1e_loop, cint_1e_drv
   public :: int1e_ovlp_cart, int1e_ovlp_sph
   public :: int1e_nuc_cart, int1e_nuc_sph
   ! int2c2e reuses the 1e output layer verbatim -- same two indices, same
   ! transform, same cache estimate.
   public :: int1e_cache_size, apply_c2s_cart_1e, apply_c2s_sph_1e, apply_c2s_dset0
   public :: int1e_ovlp_optimizer
   public :: int1e_nuc_optimizer

contains

   ! The default gout: multiply the three Cartesian factors for each component.
   subroutine cint_gout1e(gout, g, idx, envs, empty)
      real(dp), intent(inout) :: gout(0:*)
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: idx(0:*)
      type(cint_env_vars), intent(in) :: envs
      integer,  intent(in)    :: empty
      integer :: n, ix, iy, iz
      if (empty /= 0) then
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
   end subroutine cint_gout1e

   ! The nuclear-attraction gout: unlike the overlap this sums over the Rys
   ! roots.  Hand-written in the C too -- int1e_nuc's line in auto_intor.cl is
   ! commented out, so it is not generated there either.
   subroutine cint_gout1e_nuc(gout, g, idx, envs, empty)
      real(dp), intent(inout) :: gout(0:*)
      real(dp), intent(inout), target :: g(0:)
      integer,  intent(in)    :: idx(0:*)
      type(cint_env_vars), intent(in) :: envs
      integer,  intent(in)    :: empty
      integer  :: n, i, gx, gy, gz
      real(dp) :: s
      if (empty /= 0) then
         do n = 0, envs%nf - 1
            gx = idx(n*3+0); gy = idx(n*3+1); gz = idx(n*3+2)
            s = 0.0_dp
            do i = 0, envs%nrys_roots - 1
               s = s + g(gx+i) * g(gy+i) * g(gz+i)
            end do
            gout(n) = s
         end do
      else
         do n = 0, envs%nf - 1
            gx = idx(n*3+0); gy = idx(n*3+1); gz = idx(n*3+2)
            s = 0.0_dp
            do i = 0, envs%nrys_roots - 1
               s = s + g(gx+i) * g(gy+i) * g(gz+i)
            end do
            gout(n) = gout(n) + s
         end do
      end if
   end subroutine cint_gout1e_nuc

   ! One shell pair: returns .true. if anything survived screening.
   function cint_1e_loop(gctr, gctroff, envs, ws, int1e_type) result(has_value)
      real(dp), intent(inout) :: gctr(0:)
      integer,  intent(in)    :: gctroff
      type(cint_env_vars), intent(inout) :: envs
      ! TARGET so the pair table can be bound to a local pointer below
      type(cint_ws), intent(inout), target :: ws
      integer,  intent(in)    :: int1e_type
      logical :: has_value

      integer :: i_sh, j_sh, i_ctr, j_ctr, i_prim, j_prim
      logical :: use_opt
      logical :: idx_cached
      integer :: ai_o, aj_o, ci_o, cj_o, n_comp, nc
      integer :: leng, lenj, leni, len0
      integer :: og, ogout, ogctri, ogctrj, oidx, ononi, ononj, osi, osj
      integer :: olmi, olmj
      integer :: ip, jp, npair, k
      real(dp) :: expcutoff, fac1i, fac1j, expij, common_factor, rr_ij
      logical  :: allempty
      ! The C keeps three flags in one array and then ALIASES them:
      !     if (j_ctr == 1) { gctri = gctrj; iempty = jempty; }
      !     if (i_ctr == 1) { gout  = gctri; gempty = iempty; }
      ! so clearing one can clear another.  That is load-bearing -- it is how
      ! the accumulate-in-place case knows whether the buffer it shares still
      ! holds anything.  Modelled here with the same array and a pair of
      ! indices rather than by trying to reason it out into separate logicals.
      integer, parameter :: E_G = 0, E_I = 1, E_J = 2
      logical :: empty(0:2)
      integer :: ig, ii, ij
      ! No aip/ajp/cip/cjp here, though the two-electron loop hoists exactly
      ! those: measured, it saves 25,040 instructions of 506 million.  The
      ! hoist pays over a four-deep primitive nest whose innermost body is
      ! small; this nest is two deep and its body is a whole 2D recursion.
      ! The three arrays the gout kernel is called with, bound once.  The C
      ! passes three bare pointers; passing ws%d(ogout:) and ws%i(oidx:) as
      ! sections makes gfortran build a descriptor at every call, and for
      ! int1e_nuc that call happens once per nucleus per primitive pair.
      real(dp), pointer, contiguous :: gp(:), goutp(:)
      integer,  pointer, contiguous :: idxp(:)

      i_sh = envs%shls(0)
      j_sh = envs%shls(1)
      i_ctr = envs%x_ctr(0)
      j_ctr = envs%x_ctr(1)
      i_prim = bas_of(envs, NPRIM_OF, i_sh)
      j_prim = bas_of(envs, NPRIM_OF, j_sh)
      ai_o = bas_of(envs, PTR_EXP,   i_sh)
      aj_o = bas_of(envs, PTR_EXP,   j_sh)
      ci_o = bas_of(envs, PTR_COEFF, i_sh)
      cj_o = bas_of(envs, PTR_COEFF, j_sh)
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor
      expcutoff = envs%expcutoff
      has_value = .false.

      npair = i_prim * j_prim
      call ws_ensure_pd(ws, npair)
      call ws_alloc_d(ws, i_prim, olmi)
      call ws_alloc_d(ws, j_prim, olmj)
      use_opt = .false.
      if (associated(ws%opt)) use_opt = cint_opt_usable(ws%opt, 2, envs)
      if (use_opt) then
         call ws_opt_log_maxc(ws, olmi, i_sh, i_prim)
         call ws_opt_log_maxc(ws, olmj, j_sh, j_prim)
      else
      call cint_log_max_pgto_coeff(ws%d(olmi:), envs%env, ci_o, i_prim, i_ctr)
      call cint_log_max_pgto_coeff(ws%d(olmj:), envs%env, cj_o, j_prim, j_ctr)
      end if

      rr_ij = cint_square_dist(envs%ri, envs%rj)
      if (cint_set_pairdata(ws%p, envs%env, ai_o, envs%env, aj_o, envs%ri, envs%rj, &
                            ws%d(olmi:), ws%d(olmj:), envs%li_ceil, envs%lj_ceil, &
                            i_prim, j_prim, rr_ij, expcutoff, &
                            envs%env(PTR_RANGE_OMEGA))) then
         return
      end if

      call ws_alloc_i(ws, envs%nf * 3, oidx)
      ! NOT `use_opt .and. ws_opt_idx(...)`: Fortran does not
      ! guarantee short-circuit evaluation, so the compiler is free to
      ! call ws_opt_idx even when use_opt is false -- and that
      ! dereferences ws%opt, which is null precisely when use_opt is
      ! false.  gfortran short-circuits at -O2 and does not at -O0, so
      ! this segfaulted only in a debug build.
      idx_cached = .false.
      if (use_opt) idx_cached = ws_opt_idx(ws, oidx, cint_opt_idx_key(envs, 2), envs%nf*3)
      if (.not. idx_cached) then
         call cint_g1e_index_xyz(ws%i(oidx:), envs)
      end if

      call ws_alloc_i(ws, i_prim, ononi)
      call ws_alloc_i(ws, j_prim, ononj)
      call ws_alloc_i(ws, i_prim*i_ctr, osi)
      call ws_alloc_i(ws, j_prim*j_ctr, osj)
      if (use_opt) then
         call ws_opt_non0(ws, ononi, osi, i_sh, i_prim, i_ctr, 0, 0)
         call ws_opt_non0(ws, ononj, osj, j_sh, j_prim, j_ctr, 0, 0)
      else
      call cint_non0coeff_byshell(ws%i(osi:), ws%i(ononi:), envs%env, ci_o, i_prim, i_ctr)
      call cint_non0coeff_byshell(ws%i(osj:), ws%i(ononj:), envs%env, cj_o, j_prim, j_ctr)
      end if

      nc = i_ctr * j_ctr
      ! (irys,i,j,coord); +1 for nabla-r12
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenj = envs%nf * nc * n_comp
      leni = envs%nf * i_ctr * n_comp
      len0 = envs%nf * n_comp

      call ws_alloc_d(ws, leng, og)
      ! The C aliases these three onto gctr and each other whenever a
      ! contraction count is 1, so the accumulation happens in place.  Same
      ! here, with offsets instead of pointers.
      if (n_comp == 1) then
         ogctrj = gctroff
      else
         call ws_alloc_d(ws, lenj, ogctrj)
      end if
      if (j_ctr == 1) then
         ogctri = ogctrj
      else
         call ws_alloc_d(ws, leni, ogctri)
      end if
      if (i_ctr == 1) then
         ogout = ogctri
      else
         call ws_alloc_d(ws, len0, ogout)
      end if

      empty = .true.
      ig = E_G; ii = E_I; ij = E_J
      if (j_ctr == 1) ii = ij      ! iempty aliases jempty
      if (i_ctr == 1) ig = ii      ! gempty aliases iempty (hence maybe jempty)

      common_factor = envs%common_factor &
                    * cint_common_fac_sp(envs%i_l) * cint_common_fac_sp(envs%j_l)

      gp(0:)    => ws%d(og:og+leng-1)
      goutp(0:) => ws%d(ogout:ogout+len0-1)
      idxp(0:)  => ws%i(oidx:oidx+envs%nf*3-1)

      allempty = .true.
      do jp = 0, j_prim - 1
         envs%aj = envs%env(aj_o + jp)
         if (j_ctr == 1) then
            fac1j = common_factor * envs%env(cj_o + jp)
         else
            fac1j = common_factor
            empty(ii) = .true.
         end if

         do ip = 0, i_prim - 1
            k = jp * i_prim + ip
            if (ws%p(k)%cceij > expcutoff) cycle
            envs%ai = envs%env(ai_o + ip)
            expij = ws%p(k)%eij
            envs%rij = ws%p(k)%rij
            if (i_ctr == 1) then
               fac1i = fac1j * envs%env(ci_o + ip) * expij
            else
               fac1i = fac1j * expij
            end if
            envs%fac = fac1i

            call make_g1e_gout(ws%d, og, goutp, gp, idxp, envs, &
                               merge(1, 0, empty(ig)), int1e_type)

            ! primitive -> contracted over i (PRIM2CTR0 in the C)
            if (i_ctr > 1) then
               if (empty(ii)) then
                  call cint_prim_to_ctr_0(ws%d, ogctri, ogout, envs%env, ci_o+ip, &
                                          envs%nf*n_comp, i_prim, i_ctr)
               else
                  call cint_prim_to_ctr_1(ws%d, ogctri, ogout, envs%env, ci_o+ip, &
                                          envs%nf*n_comp, i_prim, i_ctr, &
                                          ws%i(ononi+ip), ws%i, osi + ip*i_ctr)
               end if
            end if
            empty(ii) = .false.
            allempty = .false.
         end do

         if (.not. empty(ii)) then
            if (j_ctr > 1) then
               if (empty(ij)) then
                  call cint_prim_to_ctr_0(ws%d, ogctrj, ogctri, envs%env, cj_o+jp, &
                                          envs%nf*i_ctr*n_comp, j_prim, j_ctr)
               else
                  call cint_prim_to_ctr_1(ws%d, ogctrj, ogctri, envs%env, cj_o+jp, &
                                          envs%nf*i_ctr*n_comp, j_prim, j_ctr, &
                                          ws%i(ononj+jp), ws%i, osj + jp*j_ctr)
               end if
            end if
            empty(ij) = .false.
         end if
      end do

      if (n_comp > 1 .and. .not. empty(ij)) then
         call cint_dmat_transpose(gctr(gctroff:), ws%d(ogctrj:), envs%nf*nc, n_comp)
      end if
      has_value = .not. allempty
   end function cint_1e_loop


   ! gout, g and idx arrive already bound -- see the note at their
   ! declarations in cint_1e_loop.  buf and og stay because cint_g1e_ovlp and
   ! cint_g1e_nuc build into the workspace by offset, not through g.
   !
   ! gout and idx are assumed-size, so they are forwarded as gout(0)/idx(0):
   ! an assumed-size array cannot be passed whole, and sequence association
   ! on the first element passes exactly the address the C would.
   subroutine make_g1e_gout(buf, og, gout, g, idx, envs, empty, int1e_type)
      real(dp), intent(inout) :: buf(0:)
      real(dp), intent(inout) :: gout(0:*), g(0:)
      integer,  intent(in)    :: og, idx(0:*), empty, int1e_type
      type(cint_env_vars), intent(in) :: envs
      integer :: dummy, ia, e
      select case (int1e_type)
      case (INT1E_OVLP)
         dummy = cint_g1e_ovlp(buf, og, envs)
         ! Through the procedure pointer, not to cint_gout1e directly: that is
         ! the whole point of the slot, and calling the default kernel here
         ! makes every integral compute a plain overlap.
         call envs%f_gout(gout(0), g, idx(0), envs, empty)
      case (INT1E_RINV)
         dummy = cint_g1e_nuc(buf, og, envs, -1)
         call envs%f_gout(gout(0), g, idx(0), envs, empty)
      case (INT1E_NUC)
         ! One pass per nucleus, accumulating; only the first clears.
         e = empty
         do ia = 0, envs%natm - 1
            dummy = cint_g1e_nuc(buf, og, envs, ia)
            call envs%f_gout(gout(0), g, idx(0), envs, e)
            e = 0
         end do
      case default
         error stop "cint_1e: unknown int1e_type"
      end select
   end subroutine make_g1e_gout

   ! Drive one shell pair all the way to the caller's buffer.
   function cint_1e_drv(out, dims, envs, ws, c2s_kind, int1e_type) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      integer,  intent(in)    :: c2s_kind, int1e_type
      logical :: has_value
      integer :: dmark, imark

      integer :: nc, n_comp, ogctr, n, nout, counts(0:3), nd, ni

      nc = envs%nf * envs%x_ctr(0) * envs%x_ctr(1)
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor

      call int1e_cache_size(envs, nd, ni)
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*n_comp, ogctr)
      ! NOT zeroed.  No C driver zeroes gctr -- CINT1e_drv, CINT3c2e_drv,
      ! CINT2c2e_drv, CINT3c1e_drv and CINT1e_grids_drv all just take the
      ! buffer and thread an empty flag through the loop, so the first
      ! contraction sets rather than accumulates.  The port threads the same
      ! flag, so the memset was pure extra writes.

      has_value = cint_1e_loop(ws%d, ogctr, envs, ws, int1e_type)

      select case (c2s_kind)
      case (C2S_SPH_1E)
         counts(0) = (envs%i_l*2+1) * envs%x_ctr(0)
         counts(1) = (envs%j_l*2+1) * envs%x_ctr(1)
      case default
         counts(0) = envs%nfi * envs%x_ctr(0)
         counts(1) = envs%nfj * envs%x_ctr(1)
      end select
      counts(2) = 1
      counts(3) = 1
      nout = dims(0) * dims(1)

      if (has_value) then
         call ws_mark(ws, dmark, imark)
         do n = 0, n_comp - 1
            call ws_rewind(ws, dmark, imark)
            select case (c2s_kind)
            case (C2S_CART_1E)
               call apply_c2s_cart_1e(out(nout*n:), ws%d(ogctr+nc*n:), dims, envs)
            case (C2S_SPH_1E)
               call apply_c2s_sph_1e(out(nout*n:), ws%d(ogctr+nc*n:), dims, envs, ws)
            case default
               error stop "cint_1e_drv: spinor c2s is D9"
            end select
         end do
      else
         do n = 0, n_comp - 1
            call apply_c2s_dset0(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_1e_drv

   ! How much scratch one shell pair needs, in doubles and in integers.
   ! The C computes the same thing in int1e_cache_size and hands the caller a
   ! single `double *`; splitting it by type costs nothing and lets the
   ! integer side be integers.
   pure subroutine int1e_cache_size(envs, nd, ni)
      type(cint_env_vars), intent(in) :: envs
      integer, intent(out) :: nd, ni
      integer :: i_prim, j_prim, i_ctr, j_ctr, nc, n_comp
      integer :: leng, lenj, leni, len0, buflen

      i_prim = bas_of(envs, NPRIM_OF, envs%shls(0))
      j_prim = bas_of(envs, NPRIM_OF, envs%shls(1))
      i_ctr = envs%x_ctr(0)
      j_ctr = envs%x_ctr(1)
      nc = i_ctr * j_ctr
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor

      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenj = envs%nf * nc * n_comp
      leni = envs%nf * i_ctr * n_comp
      len0 = envs%nf * n_comp
      ! the two c2s_sph_1e scratch blocks
      buflen = envs%nfi * (2*envs%j_l + 1)

      nd = envs%nf*nc*n_comp        &  ! gctr, allocated by the driver
         + i_prim + j_prim          &  ! log_maxci, log_maxcj
         + leng + lenj + leni + len0 &
         + 2*buflen                 &
         + 64                          ! slack, so a small miscount is not fatal
      ni = envs%nf*3                &  ! idx
         + i_prim + j_prim          &  ! non0ctri, non0ctrj
         + i_prim*i_ctr + j_prim*j_ctr &
         + 64
   end subroutine int1e_cache_size

   ! Cartesian output: reorder the contracted block into the caller's layout.
   subroutine apply_c2s_cart_1e(opij, gctr, dims, envs)
      real(dp), intent(inout) :: opij(0:)
      real(dp), intent(in)    :: gctr(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(in) :: envs
      integer :: ic, jc, i, j, ni, ofj, base, gbase
      ni = dims(0)
      ofj = ni * envs%nfj
      gbase = 0
      do jc = 0, envs%x_ctr(1) - 1
         do ic = 0, envs%x_ctr(0) - 1
            base = ofj * jc + envs%nfi * ic
            do j = 0, envs%nfj - 1
               do i = 0, envs%nfi - 1
                  opij(base + j*ni + i) = gctr(gbase + j*envs%nfi + i)
               end do
            end do
            gbase = gbase + envs%nf
         end do
      end do
   end subroutine apply_c2s_cart_1e

   ! Spherical output: transform the ket index, then the bra index, then
   ! place the block.
   subroutine apply_c2s_sph_1e(opij, gctr, dims, envs, ws)
      real(dp), intent(inout) :: opij(0:)
      real(dp), intent(in)    :: gctr(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(in) :: envs
      type(cint_ws), intent(inout) :: ws

      integer :: i_l, j_l, di, dj, ni, ofj, nfi, buflen, o1, o2
      integer :: ic, jc, i, j, base, gbase, loc1, loc2

      i_l = envs%i_l
      j_l = envs%j_l
      di = i_l*2 + 1
      dj = j_l*2 + 1
      ni = dims(0)
      ofj = ni * dj
      nfi = envs%nfi
      buflen = nfi * dj
      call ws_alloc_d(ws, buflen, o1)
      call ws_alloc_d(ws, buflen, o2)

      gbase = 0
      do jc = 0, envs%x_ctr(1) - 1
         do ic = 0, envs%x_ctr(0) - 1
            loc1 = cint_c2s_ket_sph(ws%d(o1:), gctr(gbase:), nfi, nfi, j_l)
            if (loc1 == RESULT_IN_GCART) then
               loc2 = cint_c2s_bra_sph(ws%d(o2:), dj, gctr(gbase:), i_l)
            else
               loc2 = cint_c2s_bra_sph(ws%d(o2:), dj, ws%d(o1:), i_l)
            end if
            base = ofj * jc + di * ic
            do j = 0, dj - 1
               do i = 0, di - 1
                  if (loc2 == RESULT_IN_GCART) then
                     if (loc1 == RESULT_IN_GCART) then
                        opij(base + j*ni + i) = gctr(gbase + j*di + i)
                     else
                        opij(base + j*ni + i) = ws%d(o1 + j*di + i)
                     end if
                  else
                     opij(base + j*ni + i) = ws%d(o2 + j*di + i)
                  end if
               end do
            end do
            gbase = gbase + envs%nf
         end do
      end do
   end subroutine apply_c2s_sph_1e

   pure subroutine apply_c2s_dset0(out, dims, counts)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), counts(0:)
      integer :: i, j, ni
      ni = dims(0)
      do j = 0, counts(1) - 1
         do i = 0, counts(0) - 1
            out(j*ni + i) = 0.0_dp
         end do
      end do
   end subroutine apply_c2s_dset0

   function int1e_nuc_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout1e_nuc
      has_value = cint_1e_drv(out, dims, envs, ws, C2S_CART_1E, INT1E_NUC)
   end function int1e_nuc_cart

   function int1e_nuc_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout1e_nuc
      has_value = cint_1e_drv(out, dims, envs, ws, C2S_SPH_1E, INT1E_NUC)
   end function int1e_nuc_sph

   ! The int1e_ovlp entry points that used to live here were hand-written for
   ! D5, to prove the driver.  They are now generated along with everything
   ! else -- see fortran/src/cint_intor1.f90 -- and were deleted rather than
   ! left to rot beside a generated copy of themselves.

   ! int1e_ovlp and int1e_nuc are hand-written here rather than generated,
   ! because upstream comments both out of auto_intor.cl and hand-writes them
   ! in cint1e.c.  The generator could emit them now -- the rinv path arrived
   ! with D9's two-electron operations -- but matching where the C draws the
   ! line keeps the generated catalogue a description-for-description mirror
   ! of the C's, rather than a superset that has to be explained.

   function int1e_ovlp_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout1e
      has_value = cint_1e_drv(out, dims, envs, ws, C2S_CART_1E, INT1E_OVLP)
   end function int1e_ovlp_cart

   function int1e_ovlp_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int1e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout1e
      has_value = cint_1e_drv(out, dims, envs, ws, C2S_SPH_1E, INT1E_OVLP)
   end function int1e_ovlp_sph

   ! ---- optimizers ------------------------------------------------------

   subroutine int1e_ovlp_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_all_1e_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int1e_ovlp_optimizer

   subroutine int1e_nuc_optimizer(opt, atm, natm, bas, nbas, env)
      type(cint_opt_t), intent(inout) :: opt
      integer,  intent(in) :: natm, nbas
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_all_1e_optimizer(opt, ng, atm, natm, bas, nbas, env)
   end subroutine int1e_nuc_optimizer


end module cint_1e
