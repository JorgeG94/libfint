!
! The one-electron grids spinor driver, and the five integrals on it.
!
! Ported from CINT1e_grids_spinor_drv (src/cint1e_grids.c:286) and the four
! c2s_s[fi]_1e_grids[i] transforms (src/cart2sph.c:5111-5316).
!
! Underneath is the same cint_1e_grids_loop the real drivers use.  What the
! grid index changes is only the shape of the output: the transform runs once
! per GRID_BLKSIZE block rather than once per contraction pair, because that is
! how the loop lays gctr out, and the caller's array is strided by the full
! grid count while the block is not.
!
! One wrinkle that is not present in cint_1e_spinor.  There the four spinor
! transforms differ along two axes, spin-free/spin-included and with/without a
! factor of i, and all four are reachable.  Here the C generates only two of
! them -- sf for int1e_grids, _ip, _ipvip and _ipip, si for _spvsp -- and never
! emits an entry point for the two _gridsi forms.  Both are implemented anyway:
! the choice costs one branch, and leaving half a transform out is the kind of
! gap that is only discovered by the caller who needs it.
!
module cint_1e_grids_spinor
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws, ws_ensure, ws_alloc_d, ws_mark, ws_rewind
   use cint_g1e,       only: cint_init_int1e_envvars
   use cint_1e_grids,  only: cint_init_int1e_grids_envvars, cint_gout1e_grids, &
                             cint_1e_grids_loop, GRID_BLKSIZE
   use cint_bas,       only: NPRIM_OF, KAPPA_OF
   use cint_cart2spinor
   implicit none
   private

   public :: cint_1e_grids_spinor_drv, int1e_grids_spinor

contains

   ! c2s_grids_zset0 (src/cart2sph.c:4775) is c2s_zset0 with the grid axis
   ! rotated to the front, which is where the caller's layout puts it.
   pure subroutine c2s_grids_zset0(out, dims, counts)
      complex(dp), intent(inout) :: out(0:)
      integer,     intent(in)    :: dims(0:), counts(0:)
      call c2s_zset0(out, [dims(2), dims(0), dims(1), dims(3)], &
                          [counts(2), counts(0), counts(1), counts(3)])
   end subroutine c2s_grids_zset0

   ! zcopy_grids_ij (src/cart2sph.c): merge one block's real and imaginary
   ! halves into the caller's complex array, whose grid stride is the full
   ! grid count and not the block's.
   pure subroutine zcopy_grids_ij(out, ooff, gR, gI, ng, ni, mg, mi, mj)
      complex(dp), intent(inout) :: out(0:)
      real(dp),    intent(in)    :: gR(0:), gI(0:)
      integer,     intent(in)    :: ooff, ng, ni, mg, mi, mj
      integer :: i, j, m, ob, gb
      ob = ooff
      gb = 0
      do j = 0, mj - 1
         do i = 0, mi - 1
         do m = 0, mg - 1
            out(ob + i*ng + m) = cmplx(gR(gb + i*mg + m), gI(gb + i*mg + m), dp)
         end do
         end do
         ob = ob + ng * ni
         gb = gb + mg * mi
      end do
   end subroutine zcopy_grids_ij

   ! ONE ARRAY, OFFSETS.  buf is the workspace; it holds both the driver's
   ! gctr and the four split scratch blocks, so it arrives once and the
   ! non-overlapping pieces are named by offset rather than by slice.
   subroutine apply_c2s_spinor_1e_grids(out, buf, goff, o1r, o1i, o2r, o2i, &
                                        dims, envs, spin_included, with_i)
      complex(dp), intent(inout) :: out(0:)
      real(dp),    intent(inout) :: buf(0:)
      integer,     intent(in)    :: goff, o1r, o1i, o2r, o2i, dims(0:)
      type(cint_env_vars), intent(in) :: envs
      logical,     intent(in)    :: spin_included, with_i

      integer :: ngrids, i_kp, j_kp, i_ctr, j_ctr, di, dj, ni, ng, ofj
      integer :: nfj, nf, ic, jc, goc, bgrids, bg_di, bg_nf, blk
      integer :: gx, gy, gz, g1

      ngrids = envs%nfl
      i_kp = bas_of(envs, KAPPA_OF, envs%shls(0))
      j_kp = bas_of(envs, KAPPA_OF, envs%shls(1))
      i_ctr = envs%x_ctr(0); j_ctr = envs%x_ctr(1)
      di = cint_len_spinor_kl(i_kp, envs%i_l)
      dj = cint_len_spinor_kl(j_kp, envs%j_l)
      ni = dims(0)
      ng = dims(2)
      ofj = ni * dj
      nfj = envs%nfj
      nf = envs%nf

      ! the spin-included transform reads four blocks laid end to end, each
      ! the full ngrids*nf*i_ctr*j_ctr long, in the order x, y, z, scalar
      blk = ngrids * nf * i_ctr * j_ctr
      gx = goff; gy = goff + blk; gz = goff + 2*blk; g1 = goff + 3*blk

      goc = 0
      do while (goc < ngrids)
         bgrids = min(ngrids - goc, GRID_BLKSIZE)
         bg_di = bgrids * di
         bg_nf = bgrids * nf
         do jc = 0, j_ctr - 1
         do ic = 0, i_ctr - 1
            if (spin_included) then
               call a_bra1_cart2spinor_si(buf(o1r:), buf(o1i:), buf(gx:), buf(gy:), &
                                          buf(gz:), buf(g1:), bgrids, nfj, i_kp, envs%i_l)
               gx = gx + bg_nf; gy = gy + bg_nf; gz = gz + bg_nf; g1 = g1 + bg_nf
            else
               call a_bra1_cart2spinor_sf(buf(o1r:), buf(o1i:), buf(gx:), &
                                          bgrids, nfj, i_kp, envs%i_l)
               gx = gx + bg_nf
            end if
            if (with_i) then
               call a_iket_cart2spinor(buf(o2r:), buf(o2i:), buf(o1r:), buf(o1i:), &
                                       bg_di, j_kp, envs%j_l)
            else
               call a_ket_cart2spinor(buf(o2r:), buf(o2i:), buf(o1r:), buf(o1i:), &
                                      bg_di, j_kp, envs%j_l)
            end if
            call zcopy_grids_ij(out, ng*(ofj*jc + di*ic) + goc, buf(o2r:), buf(o2i:), &
                                ng, ni, bgrids, di, dj)
         end do
         end do
         goc = goc + GRID_BLKSIZE
      end do
   end subroutine apply_c2s_spinor_1e_grids

   function cint_1e_grids_spinor_drv(out, dims, envs, ws, spin_included, with_i) &
         result(has_value)
      complex(dp), intent(inout) :: out(0:)
      integer,     intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      logical,     intent(in)    :: spin_included, with_i
      logical :: has_value

      integer :: ngrids, nf, nc, ogctr, n, nout, counts(0:3)
      integer :: nd, ni, ip, jp, leng, len0, leni, lenj, buflen, dmark, imark
      integer :: o1r, o1i, o2r, o2i, di

      ngrids = envs%nfl
      nf = envs%nf
      ! unlike the real driver, ncomp_e1 multiplies into nc rather than into
      ! the component loop: the spin-included transform consumes all four of
      ! its blocks in one call
      nc = ngrids * nf * envs%x_ctr(0) * envs%x_ctr(1) * envs%ncomp_e1

      ip = bas_of(envs, NPRIM_OF, envs%shls(0))
      jp = bas_of(envs, NPRIM_OF, envs%shls(1))
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      len0 = GRID_BLKSIZE * nf * envs%ncomp_e1 * envs%ncomp_tensor
      leni = len0 * envs%x_ctr(0)
      lenj = leni * envs%x_ctr(1)
      buflen = GRID_BLKSIZE * envs%nfi * (2*envs%j_l + 1)
      ! the four split real/imaginary blocks the spinor transform wants
      di = 4*envs%i_l + 2
      nd = nc*envs%ncomp_tensor + ip + jp + leng + len0 + leni + lenj &
         + GRID_BLKSIZE*(10 + envs%nrys_roots) + 2*buflen &
         + 4 * GRID_BLKSIZE * di * (envs%nfj + envs%nfj) + 256
      ni = nf*3 + ip + jp + ip*envs%x_ctr(0) + jp*envs%x_ctr(1) + 64
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*envs%ncomp_tensor, ogctr)
      ws%d(ogctr:ogctr+nc*envs%ncomp_tensor-1) = 0.0_dp

      has_value = cint_1e_grids_loop(ws%d, ogctr, envs, ws)

      counts(0) = cint_cgto_spinor_kl(envs, 0)
      counts(1) = cint_cgto_spinor_kl(envs, 1)
      counts(2) = ngrids
      counts(3) = 1
      nout = dims(0) * dims(1) * dims(2)

      if (has_value) then
         call ws_mark(ws, dmark, imark)
         buflen = GRID_BLKSIZE * di * (envs%nfj + envs%nfj)
         do n = 0, envs%ncomp_tensor - 1
            call ws_rewind(ws, dmark, imark)
            call ws_alloc_d(ws, buflen, o1r)
            call ws_alloc_d(ws, buflen, o1i)
            call ws_alloc_d(ws, buflen, o2r)
            call ws_alloc_d(ws, buflen, o2i)
            call apply_c2s_spinor_1e_grids(out(nout*n:), ws%d, ogctr + nc*n, &
                                           o1r, o1i, o2r, o2i, dims, envs, &
                                           spin_included, with_i)
         end do
      else
         do n = 0, envs%ncomp_tensor - 1
            call c2s_grids_zset0(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_1e_grids_spinor_drv

   ! CINTcgto_spinor for shell slot n, without pulling in cint_bas's own
   ! wrapper: the kappa and the l are already sitting in envs.
   pure integer function cint_cgto_spinor_kl(envs, n) result(r)
      type(cint_env_vars), intent(in) :: envs
      integer, intent(in) :: n
      integer :: l
      if (n == 0) then
         l = envs%i_l
      else
         l = envs%j_l
      end if
      r = cint_len_spinor_kl(bas_of(envs, KAPPA_OF, envs%shls(n)), l) &
        * envs%x_ctr(n)
   end function cint_cgto_spinor_kl

   ! ---- entry points ---------------------------------------------------

   function int1e_grids_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 0, 1]
      call cint_init_int1e_grids_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout1e_grids
      hv = cint_1e_grids_spinor_drv(out, dims, envs, ws, .false., .false.)
   end function int1e_grids_spinor

end module cint_1e_grids_spinor
