!
! Two-electron Coulomb integrals by rotated-axis McMurchie-Davidson.
!
! A second way to the same numbers int2e_cart/int2e_sph produce by Rys
! quadrature, for the low angular momenta where the rotated frame pays:
! with z along A->B and CD in the xz plane, every bra pair's product centre
! is a scalar on the z axis and every ket pair's is a point in a plane, so
! the Hermite expansions in x and y collapse to their ladder terms and the
! primitive loop of an (sp|sp) class is a handful of multiply-adds on top
! of the Boys function.  doc/ROT_AXIS_MMD.md has the derivation; the
! kernels themselves are generated from it by scripts/rotaxis_mmd.
!
! This driver takes exactly the shell-quartet interface of int2e_cart: the
! same shls/atm/bas/env, the same dims-strided output block, the same
! contraction-major component layout, L shells included.  What it does is
!
!   1. permute the quartet into the canonical class order the generated
!      kernels are written for -- shell KINDS s < p < L < d ascending within
!      each pair and between the pairs, where an L shell is one four-
!      component slot whose s and p share the primitive loops;
!   2. build the local frame and the two primitive-pair tables;
!   3. run the kernel, which returns the local-frame Cartesian block;
!   4. rotate each index back to the lab frame -- composed with the
!      Cartesian-to-spherical transform when a spherical block was asked
!      for -- and scatter into the caller's layout.
!
! It is NOT bit-identical to the Rys path and cannot be: it is a different
! algorithm.  rotaxis_check holds it to the Rys path at a tolerance that is
! documented there.
!
module cint_rotaxis_2e
   use cint_const, only: dp
   use cint_bas, only: cint_len_cart, cint_bas_is_sp, &
                       ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, &
                       PTR_COORD, ATM_SLOTS, BAS_SLOTS, NF_SP
   use cint_g1e, only: cint_common_fac_sp
   use cint_envs, only: PTR_EXPCUTOFF, EXPCUTOFF, MIN_EXPCUTOFF
   use cint_tab_cart2sph, only: g_trans_cart2sph
   use cint_workspace, only: cint_ws
   use cint_rotaxis_kernels, only: rotaxis_kernel, rotaxis_has_class, ROTAXIS_LMAX, &
                                   ROTAXIS_KIND_S, ROTAXIS_KIND_P, ROTAXIS_KIND_L, ROTAXIS_KIND_D
   implicit none
   private

   public :: int2e_rotaxis_cart, int2e_rotaxis_sph, rotaxis_supported

   ! Start of each l's block in g_trans_cart2sph (the (m, f) layout at
   ! offset + m*nf + f), as cint_cart2sph has it.
   integer, parameter :: C2S_OFFSET(0:15) = [ &
      0, 1, 10, 40, 110, 245, 476, 840, 1380, 2145, 3190, 4576, 6370, 8645, 11480, 14960 ]

   real(dp), parameter :: PI = 3.14159265358979323846_dp
   real(dp), parameter :: TWO_PI_52 = 2.0_dp * PI**2 * sqrt(PI)

   ! One shell of the quartet, as the driver sees it.  `rank` is the kind
   ! rank the kernels are ordered by; `ntype` the number of coefficient
   ! column blocks (2 for an L shell: s columns then p columns).
   type shell_t
      integer  :: l, nprim, nctr, pe, pc, rank, ntype
      logical  :: sp
      real(dp) :: r(3)
   end type shell_t

contains

   ! Can this quartet go through the rotated-axis path at all?  An L shell
   ! counts as p for the purpose.
   pure logical function rotaxis_supported(shls, bas) result(yes)
      integer, intent(in) :: shls(0:), bas(0:)
      integer :: x, l(0:3)
      do x = 0, 3
         l(x) = bas(BAS_SLOTS*shls(x) + ANG_OF)
         if (cint_bas_is_sp(shls(x), bas)) l(x) = 1
      end do
      yes = all(l <= ROTAXIS_LMAX)
   end function rotaxis_supported

   function int2e_rotaxis_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  intent(in)    :: atm(0:), bas(0:)
      real(dp), intent(in)    :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      has_value = rotaxis_drv(out, dims, shls, atm, natm, bas, nbas, env, .false.)
   end function int2e_rotaxis_cart

   function int2e_rotaxis_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  intent(in)    :: atm(0:), bas(0:)
      real(dp), intent(in)    :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      has_value = rotaxis_drv(out, dims, shls, atm, natm, bas, nbas, env, .true.)
   end function int2e_rotaxis_sph

   ! ---------------------------------------------------------------------

   function rotaxis_drv(out, dims, shls, atm, natm, bas, nbas, env, sph) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  intent(in)    :: atm(0:), bas(0:)
      real(dp), intent(in)    :: env(0:)
      logical,  intent(in)    :: sph
      logical :: has_value

      type(shell_t) :: sh(0:3)
      integer  :: x, ia, pr, ncomp_out(0:3), dim(0:3)
      integer  :: i, j, k, l
      real(dp) :: cutoff

      do x = 0, 3
         sh(x)%l     = bas(BAS_SLOTS*shls(x) + ANG_OF)
         sh(x)%nprim = bas(BAS_SLOTS*shls(x) + NPRIM_OF)
         sh(x)%nctr  = bas(BAS_SLOTS*shls(x) + NCTR_OF)
         sh(x)%pe    = bas(BAS_SLOTS*shls(x) + PTR_EXP)
         sh(x)%pc    = bas(BAS_SLOTS*shls(x) + PTR_COEFF)
         sh(x)%sp    = cint_bas_is_sp(shls(x), bas)
         ia = bas(BAS_SLOTS*shls(x) + ATOM_OF)
         pr = atm(ATM_SLOTS*ia + PTR_COORD)
         sh(x)%r = env(pr:pr+2)
         sh(x)%ntype = 1
         if (sh(x)%sp) then
            ! s and p on shared exponents: coefficient columns 0..nctr-1 are
            ! the s, nctr..2nctr-1 the p, and the output is [s p p p] per
            ! contraction, spherical or not.
            sh(x)%ntype = 2
            sh(x)%rank = ROTAXIS_KIND_L
            ncomp_out(x) = NF_SP
         else
            select case (sh(x)%l)
            case (0); sh(x)%rank = ROTAXIS_KIND_S
            case (1); sh(x)%rank = ROTAXIS_KIND_P
            case (2); sh(x)%rank = ROTAXIS_KIND_D
            case default
               error stop "cint_rotaxis_2e: angular momentum beyond the generated kernels"
            end select
            if (sph) then
               ncomp_out(x) = 2*sh(x)%l + 1
            else
               ncomp_out(x) = cint_len_cart(sh(x)%l)
            end if
         end if
         dim(x) = ncomp_out(x) * sh(x)%nctr
      end do

      if (env(PTR_EXPCUTOFF) == 0.0_dp) then
         cutoff = EXPCUTOFF
      else
         cutoff = max(MIN_EXPCUTOFF, env(PTR_EXPCUTOFF)) + 1.0_dp
      end if

      call one_block(out, dims, sh, ncomp_out, env, cutoff, sph, has_value)
      ! the scatter writes every element of the block; only a screened-out
      ! quartet needs zeroing
      if (.not. has_value) then
         do l = 0, dim(3) - 1
            do k = 0, dim(2) - 1
               do j = 0, dim(1) - 1
                  do i = 0, dim(0) - 1
                     out(i + dims(0)*(j + dims(1)*(k + dims(2)*l))) = 0.0_dp
                  end do
               end do
            end do
         end do
      end if
   end function rotaxis_drv

   ! The quartet, computed in canonical order and scattered into `out`.
   !
   ! No allocation here: the pair tables, the block and its transform
   ! scratch are automatic arrays sized by the quartet, which on a
   ! (pp|pp) of a few primitives is a few kB of stack.  The first version
   ! allocated all of them per call and that, not the kernel, was the
   ! per-quartet floor.
   subroutine one_block(out, dims, sh, ncomp_out, env, cutoff, sph, has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:)
      type(shell_t), intent(in) :: sh(0:3)
      integer,  intent(in)    :: ncomp_out(0:3)
      real(dp), intent(in)    :: env(0:), cutoff
      logical,  intent(in)    :: sph
      logical,  intent(out)   :: has_value

      integer  :: perm(0:3), m, t, code
      integer  :: nf(0:3), nout(0:3), nctr(0:3), rk(0:3)
      real(dp) :: rot(3, 3), gc(3), cl(3), rab, rcd, tmat(6, 6)
      real(dp) :: fc(0:3, 0:9)
      integer  :: nbra, nket, ncb, nck, ncomp, ncol, nb, na, npb, npk
      integer  :: ia, ib, ic, id, cca, ccb, ccc, ccd, ostr(0:3)
      integer  :: sa, sb, sc, sd, oa, ob, oc, od, xd, yc, tb, ub, uc, ud, t0
      real(dp) :: fb, fcc, fd

      ! 1. canonical order: ranks ascend within each pair, then the pairs.
      perm = [0, 1, 2, 3]
      rk = sh%rank
      if (rk(perm(0)) > rk(perm(1))) call iswap(perm(0), perm(1))
      if (rk(perm(2)) > rk(perm(3))) call iswap(perm(2), perm(3))
      if (rk(perm(0)) > rk(perm(2)) .or. &
          (rk(perm(0)) == rk(perm(2)) .and. rk(perm(1)) > rk(perm(3)))) then
         call iswap(perm(0), perm(2))
         call iswap(perm(1), perm(3))
      end if
      code = 0
      do m = 0, 3
         code = 10*code + rk(perm(m))
         nctr(m) = sh(perm(m))%nctr
         nout(m) = ncomp_out(perm(m))
         if (sh(perm(m))%sp) then
            nf(m) = NF_SP
            ! the s component carries the s normalisation, the p the p's
            fc(m, 0) = cint_common_fac_sp(0)
            fc(m, 1:3) = cint_common_fac_sp(1)
         else
            nf(m) = cint_len_cart(sh(perm(m))%l)
            fc(m, 0:nout(m)-1) = cint_common_fac_sp(sh(perm(m))%l)
         end if
      end do
      if (.not. rotaxis_has_class(code)) then
         error stop "cint_rotaxis_2e: quartet beyond the generated classes"
      end if

      ncb = nctr(0)*nctr(1)
      nck = nctr(2)*nctr(3)
      npb = sh(perm(0))%nprim * sh(perm(1))%nprim
      npk = sh(perm(2))%nprim * sh(perm(3))%nprim
      ncomp = nf(0)*nf(1)*nf(2)*nf(3)
      ncol  = ncb*nck
      call body(sh(perm(0))%ntype*sh(perm(1))%ntype*ncb, npb, &
                sh(perm(2))%ntype*sh(perm(3))%ntype*nck, npk, ncomp*ncol)

   contains

      ! Split off so the automatic arrays take their sizes from the
      ! quartet; the host's variables are all visible here.
      subroutine body(wb, npb, wk, npk, nblk)
         integer, intent(in) :: wb, npb, wk, npk, nblk
         real(dp) :: bp(5, npb), kab(wb, npb), kp(7, npk), kcd(wk, npk)
         real(dp) :: blk(nblk), work(nblk)

         ! 2. frame and pair tables
         call local_frame(sh(perm(0))%r, sh(perm(1))%r, sh(perm(2))%r, sh(perm(3))%r, &
                          rot, rab, rcd, gc, cl)
         call bra_pairs(sh(perm(0)), sh(perm(1)), env, rab, cutoff, bp, kab, nbra)
         call ket_pairs(sh(perm(2)), sh(perm(3)), env, rcd, gc, cl, cutoff, kp, kcd, nket)
         has_value = .false.
         if (nbra == 0 .or. nket == 0) return

         ! 3. the kernel: the local-frame block, (a,b,c,d) with a fastest, then
         !    the contraction columns (cA, cB, cC, cD) with cA fastest.  The
         !    2 pi^(5/2) rides in the ket weights.
         call rotaxis_kernel(code, nbra, ncb, bp, kab, nket, nck, kp, kcd, gc, cutoff, blk, has_value)
         if (.not. has_value) return

         ! 4. rotate each index back to the lab frame, composed with the
         !    spherical transform when asked for.  An s index is the identity
         !    and is skipped.
         nb = 1
         na = nf(1)*nf(2)*nf(3)*ncol
         do m = 0, 3
            if (nf(m) > 1) then
               call index_matrix(sh(perm(m))%l, sh(perm(m))%sp, rot, sph, tmat, nout(m))
               call transform_index(blk, nb, nf(m), nout(m), na, tmat, work)
            end if
            nb = nb * nout(m)
            if (m < 3) na = na / nf(m+1)
         end do

         ! 5. scatter.  Canonical position m holds original slot perm(m); the
         !    per-l normalisation goes on here because for an L shell it is
         !    per component.
         ! THE SCATTER, walked component-then-contraction so that no index
         ! needs a division.  Written the obvious way -- one flat loop per
         ! slot, splitting the index back out with mod and / -- it costs
         ! eight integer divisions per output element, and profiled inside a
         ! Fock build this routine was the largest single item in the run,
         ! larger than any integral kernel.  The strides and the
         ! normalisation products are carried down the nest instead.
         ostr(0) = 1
         ostr(1) = dims(0)
         ostr(2) = dims(0)*dims(1)
         ostr(3) = dims(0)*dims(1)*dims(2)
         sa = ostr(perm(0)); sb = ostr(perm(1))
         sc = ostr(perm(2)); sd = ostr(perm(3))
         do ccd = 0, nctr(3) - 1
            ud = nctr(2)*ccd
            do id = 0, nout(3) - 1
               od = (id + nout(3)*ccd)*sd
               xd = nout(2)*id
               fd = fc(3, id)
               do ccc = 0, nctr(2) - 1
                  uc = nctr(1)*(ccc + ud)
                  do ic = 0, nout(2) - 1
                     oc = od + (ic + nout(2)*ccc)*sc
                     yc = nout(1)*(ic + xd)
                     fcc = fd*fc(2, ic)
                     do ccb = 0, nctr(1) - 1
                        ub = nctr(0)*(ccb + uc)
                        do ib = 0, nout(1) - 1
                           ob = oc + (ib + nout(1)*ccb)*sb
                           tb = nout(0)*(ib + yc)
                           fb = fcc*fc(1, ib)
                           do cca = 0, nctr(0) - 1
                              t0 = tb + nb*(cca + ub)
                              oa = ob + nout(0)*cca*sa
                              do ia = 0, nout(0) - 1
                                 out(oa + ia*sa) = blk(t0 + ia + 1)*fb*fc(0, ia)
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end subroutine body
   end subroutine one_block

   pure subroutine iswap(a, b)
      integer, intent(inout) :: a, b
      integer :: t
      t = a; a = b; b = t
   end subroutine iswap

   ! The local frame: columns of `rot` are the local x, y, z axes in lab
   ! coordinates, so x_lab = rot * x_loc.  z is along A->B; x is the part
   ! of C->D perpendicular to z, so CD lies in the xz plane; y = z cross x.
   ! gc holds (Rs, Rc, cy); cl is C in local coordinates.
   subroutine local_frame(ra, rb, rc, rd, rot, rab, rcd, gc, cl)
      real(dp), intent(in)  :: ra(3), rb(3), rc(3), rd(3)
      real(dp), intent(out) :: rot(3, 3), rab, rcd, gc(3), cl(3)
      real(dp) :: zh(3), xh(3), yh(3), ab(3), cd(3), w(3), wn, dl(3)
      integer  :: kmin

      ab = rb - ra
      rab = sqrt(sum(ab*ab))
      if (rab > 0.0_dp) then
         zh = ab / rab
      else
         zh = [0.0_dp, 0.0_dp, 1.0_dp]
      end if
      cd = rd - rc
      rcd = sqrt(sum(cd*cd))
      w = cd - sum(cd*zh)*zh
      wn = sqrt(sum(w*w))
      if (wn > 1.0e-12_dp * max(rcd, 1.0_dp)) then
         xh = w / wn
      else
         ! CD parallel to AB (or a point): any x perpendicular to z will
         ! do.  Take the lab axis least aligned with z and project.
         kmin = minloc(abs(zh), 1)
         w = 0.0_dp; w(kmin) = 1.0_dp
         w = w - sum(w*zh)*zh
         xh = w / sqrt(sum(w*w))
      end if
      yh = [zh(2)*xh(3) - zh(3)*xh(2), zh(3)*xh(1) - zh(1)*xh(3), zh(1)*xh(2) - zh(2)*xh(1)]
      rot(:, 1) = xh; rot(:, 2) = yh; rot(:, 3) = zh

      cl = [sum(xh*(rc - ra)), sum(yh*(rc - ra)), sum(zh*(rc - ra))]
      dl = [sum(xh*(rd - ra)), sum(yh*(rd - ra)), sum(zh*(rd - ra))]
      gc(1) = dl(1) - cl(1)        ! Rs = R_CD sin(gamma)
      gc(2) = dl(3) - cl(3)        ! Rc = R_CD cos(gamma)
      gc(3) = cl(2)                ! cy
   end subroutine local_frame

   ! The bra pair table and its weights.  kab is (ntt*ncb, nbra): for the
   ! contraction pair cc = ci + nctr_a*cj and the coefficient-type pair
   ! tt = ta + ntype_a*tb, the weight c_a(ta,ci) c_b(tb,cj) exp(-mu R^2)/p
   ! sits at tt + ntt*cc + 1; a plain shell has one type and the layout is
   ! kab(cc+1, pair).
   subroutine bra_pairs(sa, sb, env, rab, cutoff, bp, kab, nbra)
      type(shell_t), intent(in) :: sa, sb
      real(dp), intent(in)  :: env(0:), rab, cutoff
      real(dp), intent(out) :: bp(5, *), kab(sa%ntype*sb%ntype*sa%nctr*sb%nctr, *)
      integer,  intent(out) :: nbra
      integer  :: ia, ib, ci, cj, ta, tb, ntt
      real(dp) :: a, b, p, eab, e

      ntt = sa%ntype * sb%ntype
      nbra = 0
      do ib = 0, sb%nprim - 1
         b = env(sb%pe + ib)
         do ia = 0, sa%nprim - 1
            a = env(sa%pe + ia)
            p = a + b
            eab = a*b/p * rab*rab
            if (eab > cutoff) cycle
            nbra = nbra + 1
            bp(1, nbra) = p
            bp(2, nbra) = 0.5_dp/p
            bp(3, nbra) = b*rab/p
            bp(4, nbra) = b*rab/p - rab
            bp(5, nbra) = eab
            e = exp(-eab)/p
            do cj = 0, sb%nctr - 1
               do ci = 0, sa%nctr - 1
                  do tb = 0, sb%ntype - 1
                     do ta = 0, sa%ntype - 1
                        kab(ta + sa%ntype*tb + ntt*(ci + sa%nctr*cj) + 1, nbra) = e &
                           * env(sa%pc + (ta*sa%nctr + ci)*sa%nprim + ia) &
                           * env(sb%pc + (tb*sb%nctr + cj)*sb%nprim + ib)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine bra_pairs

   subroutine ket_pairs(sc, sd, env, rcd, gc, cl, cutoff, kp, kcd, nket)
      type(shell_t), intent(in) :: sc, sd
      real(dp), intent(in)  :: env(0:), rcd, gc(3), cl(3), cutoff
      real(dp), intent(out) :: kp(7, *), kcd(sc%ntype*sd%ntype*sc%nctr*sd%nctr, *)
      integer,  intent(out) :: nket
      integer  :: ic, id, ck, jl, tc, td, ntt
      real(dp) :: c, d, q, ecd, e, yd

      ntt = sc%ntype * sd%ntype
      nket = 0
      do id = 0, sd%nprim - 1
         d = env(sd%pe + id)
         do ic = 0, sc%nprim - 1
            c = env(sc%pe + ic)
            q = c + d
            ecd = c*d/q * rcd*rcd
            if (ecd > cutoff) cycle
            nket = nket + 1
            yd = d/q
            kp(1, nket) = q
            kp(2, nket) = 0.5_dp/q
            kp(3, nket) = c/q
            kp(4, nket) = yd
            kp(5, nket) = cl(1) + yd*gc(1)
            kp(6, nket) = cl(3) + yd*gc(2)
            kp(7, nket) = ecd
            e = TWO_PI_52 * exp(-ecd)/q
            do jl = 0, sd%nctr - 1
               do ck = 0, sc%nctr - 1
                  do td = 0, sd%ntype - 1
                     do tc = 0, sc%ntype - 1
                        kcd(tc + sc%ntype*td + ntt*(ck + sc%nctr*jl) + 1, nket) = e &
                           * env(sc%pc + (tc*sc%nctr + ck)*sc%nprim + ic) &
                           * env(sd%pc + (td*sd%nctr + jl)*sd%nprim + id)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine ket_pairs

   ! Position of the Cartesian function x^lx y^ly z^lz in libcint's order.
   pure integer function cart_index(lx, ly, lz) result(idx)
      integer, intent(in) :: lx, ly, lz
      integer :: l, k
      l = lx + ly + lz
      idx = 0
      do k = l, lx + 1, -1
         idx = idx + (l - k + 1)
      end do
      idx = idx + (l - lx - ly)
   end function cart_index

   ! The matrix that takes the local-frame Cartesian components of one
   ! index to the lab frame -- and on to spherical when `sph` -- so that
   ! out(I) = sum_J tmat(I, J) loc(J), in the leading nout x nf of tmat.
   ! p is rot itself; L is 1 (+) rot; a lab Cartesian x^i y^j z^k of a d is
   ! the product of the corresponding rows of rot expanded as a polynomial
   ! in the local components.
   subroutine index_matrix(l, is_sp, rot, sph, tmat, nout)
      integer,  intent(in)  :: l
      logical,  intent(in)  :: is_sp
      real(dp), intent(in)  :: rot(3, 3)
      logical,  intent(in)  :: sph
      real(dp), intent(out) :: tmat(6, 6)
      integer,  intent(out) :: nout
      real(dp) :: mrot(6, 6)
      integer :: nf, nd, ix, iy, iz, row, dir, co, m, f, ff, a, b
      integer :: ex(3), fa(3), fb(3)

      if (is_sp) then
         nout = NF_SP
         tmat(1:4, 1:4) = 0.0_dp
         tmat(1, 1) = 1.0_dp
         tmat(2:4, 2:4) = rot
         return
      end if
      select case (l)
      case (1)
         nout = 3
         tmat(1:3, 1:3) = rot
         return
      case (2)
         ! x^i y^j z^k (lab) = (r_i . x)(r_j . x) summed: for the two lab
         ! factors dir1, dir2 the local monomial e_a + e_b gets rot(dir1,a)
         ! rot(dir2,b).
         nf = 6
         mrot(1:6, 1:6) = 0.0_dp
         row = 0
         do ix = 2, 0, -1
            do iy = 2 - ix, 0, -1
               iz = 2 - ix - iy
               row = row + 1
               ex = [ix, iy, iz]
               ! the two factors' directions
               fa = 0; fb = 0
               dir = 0
               do m = 1, 3
                  do f = 1, ex(m)
                     dir = dir + 1
                     if (dir == 1) fa(1) = m
                     if (dir == 2) fb(1) = m
                  end do
               end do
               do a = 1, 3
                  do b = 1, 3
                     m = cart_index(merge(1,0,a==1) + merge(1,0,b==1), &
                                    merge(1,0,a==2) + merge(1,0,b==2), &
                                    merge(1,0,a==3) + merge(1,0,b==3))
                     mrot(row, m + 1) = mrot(row, m + 1) + rot(fa(1), a) * rot(fb(1), b)
                  end do
               end do
            end do
         end do
      case default
         error stop "cint_rotaxis_2e: index_matrix beyond d"
      end select

      if (sph) then
         nd = 2*l + 1
         co = C2S_OFFSET(l)
         nout = nd
         do m = 0, nd - 1
            do ff = 1, nf
               tmat(m + 1, ff) = 0.0_dp
               do f = 0, nf - 1
                  tmat(m + 1, ff) = tmat(m + 1, ff) + g_trans_cart2sph(co + m*nf + f) * mrot(f + 1, ff)
               end do
            end do
         end do
      else
         nout = nf
         tmat(1:nf, 1:nf) = mrot(1:nf, 1:nf)
      end if
   end subroutine index_matrix

   ! blk viewed as (nb, nin, na) becomes (nb, nout, na) with the middle
   ! index transformed by tmat(nout, nin).  `work` must hold nb*nout*na.
   subroutine transform_index(blk, nb, nin, nout, na, tmat, work)
      real(dp), intent(inout) :: blk(:)
      integer,  intent(in)    :: nb, nin, nout, na
      real(dp), intent(in)    :: tmat(6, 6)
      real(dp), intent(inout) :: work(:)
      integer :: ia, io, ii, ib
      real(dp) :: acc

      do ia = 0, na - 1
         do io = 0, nout - 1
            do ib = 0, nb - 1
               acc = 0.0_dp
               do ii = 0, nin - 1
                  acc = acc + tmat(io + 1, ii + 1) * blk(ib + nb*(ii + nin*ia) + 1)
               end do
               work(ib + nb*(io + nout*ia) + 1) = acc
            end do
         end do
      end do
      blk(1:nb*nout*na) = work(1:nb*nout*na)
   end subroutine transform_index

end module cint_rotaxis_2e
