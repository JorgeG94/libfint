!
! The three-centre spinor driver, and int3c2e_spinor on top of it.
!
! Ported from CINT3c2e_spinor_drv (src/cint3c2e.c:625) and the eight
! c2s_s[fi]_3c2e1[i][_ssc] transforms (src/cart2sph.c:5997-6458).
!
! One stage, not two: with only three indices the auxiliary one is transformed
! by an ordinary spherical inner transform, and the orbital pair goes through
! the spinor bra and ket in a single pass.  The two-stage split that
! cint_2e_spinor.f90 needs is a consequence of the four-index layout, and it
! does not arise here.
!
! ssc -- "spheric-spinor-cartesian" -- leaves the auxiliary index Cartesian.
! It is the form a density-fitted relativistic code wants, because the fitting
! basis has no spin to carry, and it is the whole difference between the two
! halves of this file: the ssc transforms skip the inner spherical step and
! carry nfk where the others carry 2k+1.
!
module cint_3c2e_spinor
   use cint_const,     only: dp
   use cint_envs
   use cint_workspace, only: cint_ws, ws_ensure, ws_alloc_d, ws_mark, ws_rewind
   use cint_g2e,       only: cint_init_int3c2e_envvars, cint_init_int2c2e_envvars
   use cint_3c2e,      only: cint_3c2e_loop_nopt, cint_2c2e_loop_nopt
   use cint_1e,        only: int1e_cache_size
   use cint_1e_spinor, only: apply_c2s_spinor_1e
   use cint_2e,        only: cint_gout2e
   use cint_bas,       only: cint_cgto_spinor, cint_cgto_spheric, KAPPA_OF
   use cint_cart2sph,  only: cint_c2s_ket_sph, RESULT_IN_GCART, RESULT_IN_GSPH
   use cint_cart2spinor
   implicit none
   private

   public :: cint_3c2e_spinor_drv, int3c2e_spinor, int3c2e_spinor_ssc
   public :: cint_2c2e_spinor_drv, int2c2e_spinor

contains

   ! The two-centre spinor driver.  CINT2c2e_spinor_drv (src/cint2c2e.c:294) is
   ! the real-driver body with a complex output and the spinor transform on the
   ! way out, and it refuses -- exit(1) -- when ncomp_e1 or ncomp_e2 exceeds
   ! one.  That refusal does not bite the plain integral, which is the only
   ! two-centre spinor form libcint's own generator emits a working entry point
   ! for; the four derivative forms print "&c2s_sf_1e_spinor not implemented"
   ! and return zero instead.  Both behaviours are inherited, not invented.
   function cint_2c2e_spinor_drv(out, dims, envs, ws) result(has_value)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      logical :: has_value

      integer :: nc, n_comp, ogctr, n, nout, counts(0:3), nd, ni, dmark, imark
      logical :: empty, dummy

      if (envs%ncomp_e1 > 1 .or. envs%ncomp_e2 > 1) &
         error stop "cint_2c2e_spinor_drv not implemented"

      nc = envs%nf * envs%x_ctr(0) * envs%x_ctr(1)
      n_comp = envs%ncomp_e1 * envs%ncomp_e2 * envs%ncomp_tensor

      call int1e_cache_size(envs, nd, ni)
      ! the four split real/imaginary blocks the spinor transform allocates,
      ! which the real driver's estimate has no reason to carry
      nd = nd + 4 * envs%nfi * envs%nfj * 4 * (envs%i_l*2 + 2) + 64
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*n_comp, ogctr)
      ws%d(ogctr:ogctr+nc*n_comp-1) = 0.0_dp

      empty = .true.
      dummy = cint_2c2e_loop_nopt(ws%d, ogctr, envs, ws, empty)
      has_value = .not. empty

      counts(0) = cint_cgto_spinor(envs%shls(0), envs%bas)
      counts(1) = cint_cgto_spinor(envs%shls(1), envs%bas)
      counts(2) = 1
      counts(3) = 1
      nout = dims(0) * dims(1)

      if (has_value) then
         call ws_mark(ws, dmark, imark)
         do n = 0, n_comp - 1
            call ws_rewind(ws, dmark, imark)
            call apply_c2s_spinor_1e(out(nout*n:), ws%d, ogctr + nc*n, dims, envs, ws, &
                                     .false., .false.)
         end do
      else
         do n = 0, n_comp - 1
            call c2s_zset0(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_2c2e_spinor_drv

   function int2c2e_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int2c2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      hv = cint_2c2e_spinor_drv(out, dims, envs, ws)
   end function int2c2e_spinor

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

   ! ONE BUFFER, OFFSETS.  Same rule as the four-centre driver, and for the
   ! same reason: everything here except the caller's output lives in the
   ! workspace, so it arrives once.
   subroutine apply_c2s_spinor_3c2e1(opijk, buf, goff, o1r, o1i, o2r, o2i, ob, &
                                     dims, envs, spin_included, with_i, ssc)
      complex(dp), intent(inout) :: opijk(0:)
      real(dp),    intent(inout) :: buf(0:)
      integer,     intent(in)    :: goff, o1r, o1i, o2r, o2i, ob, dims(0:)
      type(cint_env_vars), intent(in) :: envs
      logical, intent(in) :: spin_included, with_i, ssc

      integer :: i_kp, j_kp, di, dj, dk, nfi, nfj, nfk, nf2j, nf, nfik
      integer :: ni, nj, nk, ofj, ofk, ic, jc, kc, gb, base, blk, buflen
      integer :: loc, sx, sy, sz, s1, dkk

      i_kp = bas_of(envs, KAPPA_OF, envs%shls(0))
      j_kp = bas_of(envs, KAPPA_OF, envs%shls(1))
      di = cint_len_spinor_kl(i_kp, envs%i_l)
      dj = cint_len_spinor_kl(j_kp, envs%j_l)
      dk = envs%k_l*2 + 1
      nfi = envs%nfi; nfj = envs%nfj; nfk = envs%nfk
      nf2j = nfj + nfj
      nf = envs%nf
      nfik = nfi * nfk
      ni = dims(0); nj = dims(1); nk = dims(2)
      ! ssc keeps the auxiliary index Cartesian
      if (ssc) dk = nfk
      dkk = dk
      ofj = ni * dj
      ofk = ni * nj * dkk
      buflen = nfi * dk * nfj

      blk = nf * envs%x_ctr(0) * envs%x_ctr(1) * envs%x_ctr(2)
      sx = ob; sy = ob + buflen; sz = ob + 2*buflen; s1 = ob + 3*buflen

      gb = goff
      do kc = 0, envs%x_ctr(2) - 1
      do jc = 0, envs%x_ctr(1) - 1
      do ic = 0, envs%x_ctr(0) - 1
         if (ssc) then
            ! no inner transform: the Cartesian block is already what the
            ! spinor bra wants
            if (spin_included) then
               call a_bra_cart2spinor_si(buf(o1r:), buf(o1i:), &
                                         buf(gb:), buf(gb+blk:), buf(gb+2*blk:), &
                                         buf(gb+3*blk:), nfk*nfj, i_kp, envs%i_l)
            else
               call a_bra_cart2spinor_sf(buf(o1r:), buf(o1i:), buf(gb:), &
                                         nfk*nfj, i_kp, envs%i_l)
            end if
         else if (spin_included) then
            call sph2e_inner(buf, sx, gb, envs%k_l, nfi, nfj, nfi*dk, nfik, loc)
            call sph2e_inner(buf, sy, gb+blk, envs%k_l, nfi, nfj, nfi*dk, nfik, loc)
            call sph2e_inner(buf, sz, gb+2*blk, envs%k_l, nfi, nfj, nfi*dk, nfik, loc)
            call sph2e_inner(buf, s1, gb+3*blk, envs%k_l, nfi, nfj, nfi*dk, nfik, loc)
            if (loc == RESULT_IN_GCART) then
               call a_bra_cart2spinor_si(buf(o1r:), buf(o1i:), &
                                         buf(gb:), buf(gb+blk:), buf(gb+2*blk:), &
                                         buf(gb+3*blk:), dk*nfj, i_kp, envs%i_l)
            else
               call a_bra_cart2spinor_si(buf(o1r:), buf(o1i:), &
                                         buf(sx:), buf(sy:), buf(sz:), buf(s1:), &
                                         dk*nfj, i_kp, envs%i_l)
            end if
         else
            call sph2e_inner(buf, s1, gb, envs%k_l, nfi, nfj, nfi*dk, nfik, loc)
            if (loc == RESULT_IN_GCART) then
               call a_bra_cart2spinor_sf(buf(o1r:), buf(o1i:), buf(gb:), &
                                         dk*nfj, i_kp, envs%i_l)
            else
               call a_bra_cart2spinor_sf(buf(o1r:), buf(o1i:), buf(s1:), &
                                         dk*nfj, i_kp, envs%i_l)
            end if
         end if

         if (with_i) then
            call a_iket_cart2spinor(buf(o2r:), buf(o2i:), buf(o1r:), buf(o1i:), &
                                    di*dk, j_kp, envs%j_l)
         else
            call a_ket_cart2spinor(buf(o2r:), buf(o2i:), buf(o1r:), buf(o1i:), &
                                   di*dk, j_kp, envs%j_l)
         end if
         base = ofk*kc + ofj*jc + di*ic
         ! the fourth index is a single dummy component, which is what makes
         ! the four-index reorder do a three-index one
         call zcopy_iklj(opijk, base, buf(o2r:), buf(o2i:), ni, nj, nk, di, dj, dk, 1)
         gb = gb + nf
      end do
      end do
      end do
   end subroutine apply_c2s_spinor_3c2e1

   function cint_3c2e_spinor_drv(out, dims, envs, ws, spin_included, with_i, ssc) &
         result(has_value)
      complex(dp), intent(inout) :: out(0:)
      integer,     intent(in)    :: dims(0:)
      type(cint_env_vars), intent(inout) :: envs
      type(cint_ws), intent(inout) :: ws
      logical,     intent(in)    :: spin_included, with_i, ssc
      logical :: has_value

      integer :: counts(0:3), nc, n_comp, ogctr, n, nout, nd, ni
      integer :: di, dj, dk, len1, len2, buflen, o1r, o1i, o2r, o2i, ob
      integer :: i_kp, j_kp
      logical :: empty, dummy

      i_kp = bas_of(envs, KAPPA_OF, envs%shls(0))
      j_kp = bas_of(envs, KAPPA_OF, envs%shls(1))
      di = cint_len_spinor_kl(i_kp, envs%i_l)
      dj = cint_len_spinor_kl(j_kp, envs%j_l)
      if (ssc) then
         dk = envs%nfk
      else
         dk = envs%k_l*2 + 1
      end if

      counts(0) = cint_cgto_spinor(envs%shls(0), envs%bas)
      counts(1) = cint_cgto_spinor(envs%shls(1), envs%bas)
      counts(2) = dk * envs%x_ctr(2)
      counts(3) = 1

      nc = envs%nf * envs%x_ctr(0) * envs%x_ctr(1) * envs%x_ctr(2)
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor

      len1 = di * dk * 2*envs%nfj
      len2 = di * dk * dj
      buflen = envs%nfi * dk * envs%nfj

      call int3c2e_spinor_cache_size(envs, nd, ni)
      nd = nd + 2*max(len1, len2) + 2*max(len1, len2) + 4*buflen + 256
      call ws_ensure(ws, nd, ni)
      call ws_alloc_d(ws, nc*n_comp, ogctr)
      ws%d(ogctr:ogctr+nc*n_comp-1) = 0.0_dp

      empty = .true.
      dummy = cint_3c2e_loop_nopt(ws%d, ogctr, envs, ws, empty)
      has_value = .not. empty

      nout = dims(0) * dims(1) * dims(2)
      if (has_value) then
         call ws_alloc_d(ws, max(len1, len2), o1r)
         call ws_alloc_d(ws, max(len1, len2), o1i)
         call ws_alloc_d(ws, max(len1, len2), o2r)
         call ws_alloc_d(ws, max(len1, len2), o2i)
         call ws_alloc_d(ws, 4*buflen, ob)
         do n = 0, envs%ncomp_e2 * envs%ncomp_tensor - 1
            call apply_c2s_spinor_3c2e1(out(nout*n:), ws%d, &
                                        ogctr + nc*envs%ncomp_e1*n, &
                                        o1r, o1i, o2r, o2i, ob, dims, envs, &
                                        spin_included, with_i, ssc)
         end do
      else
         do n = 0, envs%ncomp_e2 * envs%ncomp_tensor - 1
            call c2s_zset0(out(nout*n:), dims, counts)
         end do
      end if
   end function cint_3c2e_spinor_drv

   pure subroutine int3c2e_spinor_cache_size(envs, nd, ni)
      type(cint_env_vars), intent(in) :: envs
      integer, intent(out) :: nd, ni
      integer :: ip, jp, kp, ic, jc, kc, nc, n_comp, nf
      integer :: leng, lenk, lenj, leni, len0
      ip = bas_of(envs, 2, envs%shls(0))
      jp = bas_of(envs, 2, envs%shls(1))
      kp = bas_of(envs, 2, envs%shls(2))
      ic = envs%x_ctr(0); jc = envs%x_ctr(1); kc = envs%x_ctr(2)
      nf = envs%nf
      nc = ic*jc*kc
      n_comp = envs%ncomp_e1 * envs%ncomp_tensor
      leng = envs%g_size * 3 * (2**envs%gbits + 1)
      lenk = nf * nc * n_comp
      lenj = nf * ic*jc * n_comp
      leni = nf * ic * n_comp
      len0 = nf * n_comp
      nd = nf*nc*n_comp + ip+jp + leng + lenk + lenj + leni + len0 + 64
      ni = nf*3 + ip+jp+kp + ip*ic + jp*jc + kp*kc + 64
   end subroutine int3c2e_spinor_cache_size

   ! ---- entry points -------------------------------------------------

   function int3c2e_spinor(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int3c2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      hv = cint_3c2e_spinor_drv(out, dims, envs, ws, .false., .false., .false.)
   end function int3c2e_spinor

   function int3c2e_spinor_ssc(out, dims, shls, atm, natm, bas, nbas, env, ws) result(hv)
      complex(dp), intent(inout) :: out(0:)
      integer,  intent(in) :: dims(0:), shls(0:), natm, nbas
      integer,  target :: atm(0:), bas(0:)
      real(dp), target :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: hv
      type(cint_env_vars) :: envs
      integer, parameter :: ng(0:7) = [0, 0, 0, 0, 0, 1, 1, 1]
      call cint_init_int3c2e_envvars(envs, ng, shls, atm, natm, bas, nbas, env)
      envs%f_gout => cint_gout2e
      hv = cint_3c2e_spinor_drv(out, dims, envs, ws, .false., .false., .true.)
   end function int3c2e_spinor_ssc

end module cint_3c2e_spinor
