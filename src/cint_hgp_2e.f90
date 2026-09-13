!
! Two-electron Coulomb integrals by Obara-Saika with Head-Gordon-Pople's
! contract-then-transfer split.
!
! A third route to the numbers int2e_cart/int2e_sph produce.  The vertical
! recurrence builds [e0|f0] -- every unit of angular momentum sitting on A
! and C -- inside the primitive loop; those are contracted; and only then
! do the horizontal recurrences move momentum onto B and D, which they can
! do once per contracted quartet because they depend on A-B and C-D and not
! on the exponents.  That is the whole idea, and it is worth most where the
! contraction is deep and the momentum high, which is where the rotated-axis
! kernels stop paying.
!
! The recurrences are generated, not written: scripts/hgp emits one module
! per class into src/hgp/, and scripts/hgp/check_numeric.py holds them to a
! McMurchie-Davidson reference in Python before any Fortran exists.  Nothing
! is expanded into a polynomial -- each intermediate is computed once from
! earlier ones -- which is the difference between this and the rotated-axis
! d classes that lost to Rys on table lookups.
!
! Same shell-quartet interface as int2e_cart, L shells included.  It is a
! different algorithm from Rys and so cannot be bit-identical; hgp_check
! holds it to the Rys path at a documented tolerance.
!
module cint_hgp_2e
   use cint_const, only: dp
   use cint_bas, only: cint_len_cart, cint_bas_is_sp, &
                       ATOM_OF, ANG_OF, NPRIM_OF, NCTR_OF, PTR_EXP, PTR_COEFF, &
                       PTR_COORD, ATM_SLOTS, BAS_SLOTS, NF_SP
   use cint_g1e, only: cint_common_fac_sp
   use cint_envs, only: PTR_EXPCUTOFF, EXPCUTOFF, MIN_EXPCUTOFF
   use cint_tab_cart2sph, only: g_trans_cart2sph
   use cint_workspace, only: cint_ws
   use cint_hgp_kernels, only: hgp_kernel, hgp_has_class, HGP_LMAX, &
                               HGP_KIND_S, HGP_KIND_P, HGP_KIND_L, HGP_KIND_D
   implicit none
   private

   public :: int2e_hgp_cart, int2e_hgp_sph, hgp_supported

   ! Start of each l's block in g_trans_cart2sph, as cint_cart2sph has it.
   integer, parameter :: C2S_OFFSET(0:15) = [ &
      0, 1, 10, 40, 110, 245, 476, 840, 1380, 2145, 3190, 4576, 6370, 8645, 11480, 14960 ]

   type shell_t
      integer  :: l, nprim, nctr, pe, pc, rank, ntype
      logical  :: sp
      real(dp) :: r(3)
   end type shell_t

contains

   ! Every shell of the quartet must be s, p, L or d.
   pure logical function hgp_supported(shls, bas) result(yes)
      integer, intent(in) :: shls(0:), bas(0:)
      integer :: x, l(0:3)
      do x = 0, 3
         l(x) = bas(BAS_SLOTS*shls(x) + ANG_OF)
         if (cint_bas_is_sp(shls(x), bas)) l(x) = 1
      end do
      yes = all(l <= HGP_LMAX)
   end function hgp_supported

   function int2e_hgp_cart(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  intent(in)    :: atm(0:), bas(0:)
      real(dp), intent(in)    :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      has_value = hgp_drv(out, dims, shls, atm, natm, bas, nbas, env, .false.)
   end function int2e_hgp_cart

   function int2e_hgp_sph(out, dims, shls, atm, natm, bas, nbas, env, ws) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  intent(in)    :: atm(0:), bas(0:)
      real(dp), intent(in)    :: env(0:)
      type(cint_ws), intent(inout) :: ws
      logical :: has_value
      has_value = hgp_drv(out, dims, shls, atm, natm, bas, nbas, env, .true.)
   end function int2e_hgp_sph

   ! ---------------------------------------------------------------------

   function hgp_drv(out, dims, shls, atm, natm, bas, nbas, env, sph) result(has_value)
      real(dp), intent(inout) :: out(0:)
      integer,  intent(in)    :: dims(0:), shls(0:), natm, nbas
      integer,  intent(in)    :: atm(0:), bas(0:)
      real(dp), intent(in)    :: env(0:)
      logical,  intent(in)    :: sph
      logical :: has_value

      type(shell_t) :: sh(0:3)
      integer  :: x, ia, pr, ncomp_out(0:3), dim(0:3), i, j, k, l
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
            sh(x)%ntype = 2
            sh(x)%rank = HGP_KIND_L
            ncomp_out(x) = NF_SP
         else
            select case (sh(x)%l)
            case (0); sh(x)%rank = HGP_KIND_S
            case (1); sh(x)%rank = HGP_KIND_P
            case (2); sh(x)%rank = HGP_KIND_D
            case default
               error stop "cint_hgp_2e: angular momentum beyond the generated kernels"
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
   end function hgp_drv

   ! The quartet in canonical order -- the higher kind first in each pair,
   ! which is the direction the transfers run -- then scattered back.
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
      real(dp) :: ab(3), cd(3), tmat(6, 6), fc(0:3, 0:9)
      integer  :: nbra, nket, ncb, nck, ncomp, ncol, nb, na, npb, npk
      integer  :: ia, ib, ic, id, cca, ccb, ccc, ccd, ostr(0:3)
      integer  :: sa, sb, sc, sd, oa, ob, oc, od, xd, yc, tb, ub, uc, ud, t0
      real(dp) :: fb, fcc, fd

      perm = [0, 1, 2, 3]
      rk = sh%rank
      if (rk(perm(0)) < rk(perm(1))) call iswap(perm(0), perm(1))
      if (rk(perm(2)) < rk(perm(3))) call iswap(perm(2), perm(3))
      if (rk(perm(0)) < rk(perm(2)) .or. &
          (rk(perm(0)) == rk(perm(2)) .and. rk(perm(1)) < rk(perm(3)))) then
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
            fc(m, 0) = cint_common_fac_sp(0)
            fc(m, 1:3) = cint_common_fac_sp(1)
         else
            nf(m) = cint_len_cart(sh(perm(m))%l)
            fc(m, 0:nout(m)-1) = cint_common_fac_sp(sh(perm(m))%l)
         end if
      end do
      if (.not. hgp_has_class(code)) then
         error stop "cint_hgp_2e: quartet beyond the generated classes"
      end if

      ab = sh(perm(0))%r - sh(perm(1))%r
      cd = sh(perm(2))%r - sh(perm(3))%r
      ncb = nctr(0)*nctr(1)
      nck = nctr(2)*nctr(3)
      npb = sh(perm(0))%nprim * sh(perm(1))%nprim
      npk = sh(perm(2))%nprim * sh(perm(3))%nprim
      ncomp = nf(0)*nf(1)*nf(2)*nf(3)
      ncol  = ncb*nck
      call body(sh(perm(0))%ntype*sh(perm(1))%ntype*ncb, npb, &
                sh(perm(2))%ntype*sh(perm(3))%ntype*nck, npk, ncomp*ncol)

   contains

      ! Split off so the automatic arrays take their sizes from the quartet.
      subroutine body(wb, npb, wk, npk, nblk)
         integer, intent(in) :: wb, npb, wk, npk, nblk
         real(dp) :: bp(10, npb), kab(wb, npb), kp(10, npk), kcd(wk, npk)
         real(dp) :: blk(nblk), work(nblk)

         call bra_pairs(sh(perm(0)), sh(perm(1)), env, cutoff, bp, kab, nbra)
         call bra_pairs(sh(perm(2)), sh(perm(3)), env, cutoff, kp, kcd, nket)
         has_value = .false.
         if (nbra == 0 .or. nket == 0) return

         call hgp_kernel(code, nbra, ncb, bp, kab, nket, nck, kp, kcd, &
                         ab, cd, cutoff, blk, has_value)
         if (.not. has_value) return

         ! Cartesian to spherical, per index; s, p and L need none.
         nb = 1
         na = nf(1)*nf(2)*nf(3)*ncol
         do m = 0, 3
            if (sph .and. nf(m) > 3 .and. .not. sh(perm(m))%sp) then
               call c2s_matrix(sh(perm(m))%l, tmat)
               call transform_index(blk, nb, nf(m), nout(m), na, tmat, work)
            end if
            nb = nb * nout(m)
            if (m < 3) na = na / nf(m+1)
         end do

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

   ! One pair's primitive table.  Used for the ket as well: Q, QC and the
   ! ket weight have exactly the shape P, PA and the bra weight have.
   subroutine bra_pairs(sa, sb, env, cutoff, bp, kab, nbra)
      type(shell_t), intent(in) :: sa, sb
      real(dp), intent(in)  :: env(0:), cutoff
      real(dp), intent(out) :: bp(10, *), kab(sa%ntype*sb%ntype*sa%nctr*sb%nctr, *)
      integer,  intent(out) :: nbra
      integer  :: ia, ib, ci, cj, ta, tb, ntt
      real(dp) :: a, b, p, eab, rab2, ba(3)

      ntt = sa%ntype * sb%ntype
      ba = sb%r - sa%r
      rab2 = sum(ba*ba)
      nbra = 0
      do ib = 0, sb%nprim - 1
         b = env(sb%pe + ib)
         do ia = 0, sa%nprim - 1
            a = env(sa%pe + ia)
            p = a + b
            eab = a*b/p * rab2
            if (eab > cutoff) cycle
            nbra = nbra + 1
            bp(1, nbra) = p
            bp(2, nbra) = 0.5_dp/p
            ! P - A = b (B - A)/p, and P itself for the P-Q difference
            bp(3:5, nbra) = b*ba/p
            bp(6:8, nbra) = sa%r + b*ba/p
            bp(9, nbra) = exp(-eab)/p
            bp(10, nbra) = eab
            do cj = 0, sb%nctr - 1
               do ci = 0, sa%nctr - 1
                  do tb = 0, sb%ntype - 1
                     do ta = 0, sa%ntype - 1
                        kab(ta + sa%ntype*tb + ntt*(ci + sa%nctr*cj) + 1, nbra) = &
                             env(sa%pc + (ta*sa%nctr + ci)*sa%nprim + ia) &
                           * env(sb%pc + (tb*sb%nctr + cj)*sb%nprim + ib)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end subroutine bra_pairs

   ! The (2l+1) x nf Cartesian-to-spherical block, in the leading corner.
   subroutine c2s_matrix(l, tmat)
      integer,  intent(in)  :: l
      real(dp), intent(out) :: tmat(6, 6)
      integer :: nf, nd, co, m, f
      nf = cint_len_cart(l)
      nd = 2*l + 1
      co = C2S_OFFSET(l)
      do m = 0, nd - 1
         do f = 0, nf - 1
            tmat(m + 1, f + 1) = g_trans_cart2sph(co + m*nf + f)
         end do
      end do
   end subroutine c2s_matrix

   ! blk viewed as (nb, nin, na) becomes (nb, nout, na).
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

end module cint_hgp_2e
