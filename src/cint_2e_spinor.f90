!
! The two-electron spinor driver, and int2e_spinor on top of it.
!
! Ported from CINT2e_spinor_drv (src/cint2e.c:871) and the eight
! c2s_s[fi]_2e[12][i] transforms (src/cart2sph.c:5385-5820).
!
! The two-electron transform runs in two stages, and the reason is the tensor's
! shape rather than tidiness: the Cartesian block is laid out (i,k,l,j), so the
! i and j indices are the outermost and innermost and can be transformed
! together in one pass, while k and l have to wait for a second pass that
! reorders the result into the caller's (i,j,k,l).  Stage one writes a split
! real/imaginary intermediate; stage two reads it, finishes the transform, and
! merges into the complex output.
!
! That intermediate is where the "complex at the edge, split inside" rule of
! cint_cart2spinor earns itself twice over: it is not the caller's array, it is
! never seen from outside this module, and keeping it split is what lets the
! second stage read the four sigma components as four adjacent real blocks the
! way the C does.
!
module cint_2e_spinor
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws, ws_ensure, ws_alloc_d
   use cint_g2e,       only: cint_init_int2e_envvars
   use cint_2e,        only: cint_2e_loop_nopt, cint_gout2e, int2e_cache_size
   use cint_bas,       only: cint_cgto_spinor, KAPPA_OF
   use cint_cart2spinor
   implicit none
   private

   public :: cint_2e_spinor_drv, int2e_spinor

contains

   ! Stage one: transform i and j, leaving k and l Cartesian.  Output is the
   ! split intermediate, one (R block, I block) pair per contraction.
   !
   ! ONE BUFFER, FOUR OFFSETS.  The output, the input and the two scratch
   ! blocks all live in the workspace, so they arrive as one array and four
   ! offsets rather than as separate dummies.  Passing ws%d as an intent(in)
   ! argument and ws -- which contains it -- as an intent(inout) one is the
   ! same illegal aliasing that cost this port an afternoon in D6 and again in
   ! D7, and it is exactly what a narrow test stride hides: it produced the
   ! right answer for every quartet the first sweep happened to visit.
   subroutine apply_c2s_spinor_2e1(buf, ooff, goff, o1r, o1i, envs, spin_included, with_i)
      real(dp), intent(inout) :: buf(0:)
      integer,  intent(in)    :: ooff, goff, o1r, o1i
      type(cint_env_vars), intent(in) :: envs
      logical,  intent(in) :: spin_included, with_i

      integer :: i_kp, j_kp, di, dj, nfj, nfk, nfl, nf2j, nf, no, d_i, d_j
      integer :: nctr, i, gb, ob, blk

      i_kp = bas_of(envs, KAPPA_OF, envs%shls(0))
      j_kp = bas_of(envs, KAPPA_OF, envs%shls(1))
      di = cint_len_spinor_kl(i_kp, envs%i_l)
      dj = cint_len_spinor_kl(j_kp, envs%j_l)
      nfj = envs%nfj; nfk = envs%nfk; nfl = envs%nfl
      nf2j = nfj + nfj
      nf = envs%nf
      no  = di * nfk * nfl * dj
      d_i = di * nfk * nfl
      d_j = nfk * nfl * nfj
      nctr = envs%x_ctr(0) * envs%x_ctr(1) * envs%x_ctr(2) * envs%x_ctr(3)

      blk = nf * nctr          ! stride between the four sigma components
      gb = goff; ob = ooff
      do i = 0, nctr - 1
         if (spin_included) then
            call a_bra_cart2spinor_si(buf(o1r:), buf(o1i:), &
                                      buf(gb:), buf(gb+blk:), buf(gb+2*blk:), &
                                      buf(gb+3*blk:), d_j, i_kp, envs%i_l)
         else
            call a_bra_cart2spinor_sf(buf(o1r:), buf(o1i:), buf(gb:), &
                                      d_j, i_kp, envs%i_l)
         end if
         if (with_i) then
            call a_iket_cart2spinor(buf(ob:), buf(ob+no:), buf(o1r:), buf(o1i:), &
                                    d_i, j_kp, envs%j_l)
         else
            call a_ket_cart2spinor(buf(ob:), buf(ob+no:), buf(o1r:), buf(o1i:), &
                                   d_i, j_kp, envs%j_l)
         end if
         gb = gb + nf
         ob = ob + no * 2
      end do
   end subroutine apply_c2s_spinor_2e1

   ! Stage two: transform k and l, and place the block in (i,j,k,l) order.
   ! One buffer here too: fijkl is the caller's own array, but the
   ! intermediate and the scratch are both the workspace.
   subroutine apply_c2s_spinor_2e2(fijkl, buf, ooff, o1r, o1i, o2r, o2i, dims, envs, &
                                   spin_included, with_i)
      complex(dp), intent(inout) :: fijkl(0:)
      real(dp),    intent(inout) :: buf(0:)
      integer,     intent(in)    :: ooff, o1r, o1i, o2r, o2i, dims(0:)
      type(cint_env_vars), intent(in) :: envs
      logical, intent(in) :: spin_included, with_i

      integer :: i_kp, j_kp, k_kp, l_kp, di, dj, dk, dl
      integer :: ni, nj, nk, nfk, nfl, nf2l, nop, ofj, ofk, ofl
      integer :: ic, jc, kc, lc, ob, base, blk

      i_kp = bas_of(envs, KAPPA_OF, envs%shls(0))
      j_kp = bas_of(envs, KAPPA_OF, envs%shls(1))
      k_kp = bas_of(envs, KAPPA_OF, envs%shls(2))
      l_kp = bas_of(envs, KAPPA_OF, envs%shls(3))
      di = cint_len_spinor_kl(i_kp, envs%i_l)
      dj = cint_len_spinor_kl(j_kp, envs%j_l)
      dk = cint_len_spinor_kl(k_kp, envs%k_l)
      dl = cint_len_spinor_kl(l_kp, envs%l_l)
      ni = dims(0); nj = dims(1); nk = dims(2)
      nfk = envs%nfk; nfl = envs%nfl
      nf2l = nfl + nfl
      nop = di * nfk * nfl * dj
      ofj = ni * dj
      ofk = ni * nj * dk
      ofl = ni * nj * nk * dl

      blk = nop * 2 * envs%x_ctr(0) * envs%x_ctr(1) * envs%x_ctr(2) * envs%x_ctr(3)
      ob = ooff
      do lc = 0, envs%x_ctr(3) - 1
      do kc = 0, envs%x_ctr(2) - 1
      do jc = 0, envs%x_ctr(1) - 1
      do ic = 0, envs%x_ctr(0) - 1
         if (spin_included) then
            call a_bra1_cart2spinor_zi(buf(o1r:), buf(o1i:), &
                                       buf(ob:), buf(ob+blk:), buf(ob+2*blk:), &
                                       buf(ob+3*blk:), di, nfl*dj, k_kp, envs%k_l)
         else
            call a_bra1_cart2spinor_zf(buf(o1r:), buf(o1i:), buf(ob:), &
                                       di, nfl*dj, k_kp, envs%k_l)
         end if
         if (with_i) then
            call a_iket1_cart2spinor(buf(o2r:), buf(o2i:), buf(o1r:), buf(o1i:), &
                                     di*dk, dj, l_kp, envs%l_l)
         else
            call a_ket1_cart2spinor(buf(o2r:), buf(o2i:), buf(o1r:), buf(o1i:), &
                                    di*dk, dj, l_kp, envs%l_l)
         end if
         base = ofl*lc + ofk*kc + ofj*jc + di*ic
         call zcopy_iklj(fijkl, base, buf(o2r:), buf(o2i:), ni, nj, nk, di, dj, dk, dl)
         ob = ob + nop * 2
      end do
      end do
      end do
      end do
   end subroutine apply_c2s_spinor_2e2

   ! The two stages are chosen independently, because the C picks them
   ! independently: int2e_spsp1_spinor is c2s_si_2e1 followed by c2s_sf_2e2,
   ! the sigma-dot-p pair sitting on electron one alone.
   function cint_2e_spinor_drv(out, dims, envs, ws, si1, wi1, si2, wi2) result(has_value)
      complex(dp), intent(inout) :: out(0:)
      integer,     intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      logical,     intent(in)    :: si1, wi1, si2, wi2
      logical :: has_value

      integer :: counts(0:3), nf, nc, n_comp, n1, ogctr, oopij, n, m, nout, nd, ni
      integer :: di, dj, dk, dl, o1r, o1i, o2r, o2i, len1, len2
      logical :: empty, dummy

      counts(0) = cint_cgto_spinor(envs%shls(0), envs%bas)
      counts(1) = cint_cgto_spinor(envs%shls(1), envs%bas)
      counts(2) = cint_cgto_spinor(envs%shls(2), envs%bas)
      counts(3) = cint_cgto_spinor(envs%shls(3), envs%bas)
      nf = envs%nf
      nc = nf * envs%x_ctr(0) * envs%x_ctr(1) * envs%x_ctr(2) * envs%x_ctr(3)
      n_comp = envs%ncomp_e1 * envs%ncomp_e2 * envs%ncomp_tensor
      ! the split intermediate: i and j already spinor, k and l still Cartesian
      n1 = counts(0) * envs%nfk * envs%x_ctr(2) * envs%nfl * envs%x_ctr(3) * counts(1)

      di = counts(0) / max(envs%x_ctr(0), 1)
      dj = counts(1) / max(envs%x_ctr(1), 1)
      dk = counts(2) / max(envs%x_ctr(2), 1)
      dl = counts(3) / max(envs%x_ctr(3), 1)

      call int2e_cache_size(envs, nd, ni)
      ! opij, plus the two-stage scratch the spinor transforms want on top of
      ! whatever the Cartesian loop needed
      nd = nd + 2*n1*envs%ncomp_e2 &
              + 4 * max(di * envs%nfk * envs%nfl * 2*envs%nfj, &
                        di * dk * 2*envs%nfl * dj) + 256
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*n_comp, ogctr)
      ws%d(ogctr:ogctr+nc*n_comp-1) = 0.0_dp

      empty = .true.
      dummy = cint_2e_loop_nopt(ws%d, ogctr, envs, ws, empty)
      has_value = .not. empty

      nout = dims(0) * dims(1) * dims(2) * dims(3)
      if (has_value) then
         ! Everything the two stages need is allocated here, once, so the
         ! stages themselves take offsets into one buffer and no second path
         ! to it.
         call ws_alloc_d(ws, 2*n1*envs%ncomp_e2, oopij)
         ! Stage one wants two blocks of len1, stage two four of len2, and
         ! stage two reuses the first pair once stage one is done.  Four of
         ! the larger size, rather than a bespoke allocation per stage:
         ! getting the two sizes crossed writes one block into the next, and
         ! looks exactly like a transform bug.
         len1 = di * envs%nfk * envs%nfl * 2*envs%nfj
         len2 = di * dk * 2*envs%nfl * dj
         len1 = max(len1, len2)
         call ws_alloc_d(ws, len1, o1r)
         call ws_alloc_d(ws, len1, o1i)
         call ws_alloc_d(ws, len1, o2r)
         call ws_alloc_d(ws, len1, o2i)
         do n = 0, envs%ncomp_tensor - 1
            do m = 0, envs%ncomp_e2 - 1
               call apply_c2s_spinor_2e1(ws%d, oopij + 2*n1*m, &
                                         ogctr + nc*envs%ncomp_e1*(n*envs%ncomp_e2 + m), &
                                         o1r, o1i, envs, si1, wi1)
            end do
            ! stage two wants four blocks of its own size; o1r/o1i are free
            ! again by now, and o2r/o2i were sized for the larger of the two
            call apply_c2s_spinor_2e2(out(nout*n:), ws%d, oopij, &
                                      o1r, o1i, o2r, o2i, dims, envs, si2, wi2)
         end do
      else
         do n = 0, envs%ncomp_tensor - 1
            call c2s_zset0(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_2e_spinor_drv

   ! ---- entry points -------------------------------------------------

   function int2e_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  target        :: atm(0:), bas(0:)
      real(dp), target        :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      has_value = cint_2e_spinor_drv(out, dims, envs, ws, .false., .false., .false., .false.)
   end function int2e_spinor

end module cint_2e_spinor
