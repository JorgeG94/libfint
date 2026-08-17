!
! Cartesian-to-spinor transforms.
!
! Ported from the second half of src/cart2sph.c (3920-4490, 6489-6960): the
! `a_*_cart2spinor_*` kernels the c2s drivers call, and the `CINTc2s_*_spinor*`
! entry points the generated relativistic integrals call directly.
!
! The coefficients are g_trans_cart2jR and g_trans_cart2jI, ported and verified
! bitwise in D2.  Two sets per l, selected by the sign of kappa: `gt_l` is
! j = l + 1/2 (kappa < 0) and `lt_l` is j = l - 1/2 (kappa > 0).  Both stop at
! l = 12, which is where the C's table stops too.
!
! COMPLEX AT THE EDGE, SPLIT INSIDE.  The public routines take and return
! complex(dp), because that is the honest type for a spinor integral and the
! reason the C does not use it -- C's double complex optimises badly -- does
! not apply here (PORT_TO_FORTRAN.md 4.5).  Inside each kernel the real and
! imaginary parts are still accumulated separately, in the C's order, term for
! term.  That is not timidity: these sums are not associative in IEEE, and
! writing `sum(conjg(c) * g)` instead would be a different number in the last
! bit.  The C's own bra transforms carry the conjugation as sign flips on the
! imaginary terms, and reproducing them literally is what makes elementwise
! agreement checkable rather than approximate.
!
module cint_cart2spinor
   use cint_const, only: dp
   use cint_tab_cart2sph, only: g_trans_cart2jR, g_trans_cart2jI, len_cart
   implicit none
   private

   public :: cint_len_spinor_kl
   public :: c2s_bra_spinor_e1sf, c2s_bra_spinor_sf, c2s_bra_spinor
   public :: c2s_bra_spinor_si, c2s_ket_spinor, c2s_iket_spinor
   public :: c2s_ket_spinor_sf1, c2s_iket_spinor_sf1
   public :: c2s_ket_spinor_si1, c2s_iket_spinor_si1
   public :: a_bra_cart2spinor_si, a_bra_cart2spinor_sf
   public :: a_ket_cart2spinor_si, a_ket_cart2spinor_sf
   public :: a_ket_cart2spinor, a_iket_cart2spinor
   public :: a_ket1_cart2spinor, a_iket1_cart2spinor
   public :: a_bra1_cart2spinor_si, a_bra1_cart2spinor_sf
   public :: a_bra1_cart2spinor_zi, a_bra1_cart2spinor_zf
   public :: zcopy_iklj
   public :: c2s_zset0

   ! Where each l's block starts in g_trans_cart2j{R,I}.  Two ladders, because
   ! the two kappa signs have different block sizes: lt_l has 2*(2l) rows of
   ! nf, gt_l has 2*(2l+2).  These are the offsets the C spells out in the
   ! g_c2s table literal.
   integer, parameter :: OFF_LT(0:12) = &
      [0, 4, 40, 160, 440, 980, 1904, 3360, 5520, 8580, 12760, 18304, 25480]
   integer, parameter :: OFF_GT(0:12) = &
      [0, 16, 88, 280, 680, 1400, 2576, 4368, 6960, 10560, 15400, 21736, 29848]
   integer, parameter :: SPINOR_LMAX = 12

contains

   ! Components of a spinor shell.  cint_bas has the same function keyed on a
   ! shell id; this one takes kappa and l directly, which is what the
   ! transforms have to hand.
   pure function cint_len_spinor_kl(kappa, l) result(n)
      integer, intent(in) :: kappa, l
      integer :: n
      if (kappa == 0) then
         n = 4*l + 2
      else if (kappa < 0) then
         n = 2*l + 2
      else
         n = 2*l
      end if
   end function cint_len_spinor_kl

   ! Start of this (kappa, l)'s coefficient block.  Split out because every
   ! kernel below opens with it, and because it is the one place the l = 12
   ! ceiling can be stated once.
   pure subroutine coeff_offset(kappa, l, off)
      integer, intent(in)  :: kappa, l
      integer, intent(out) :: off
      if (l > SPINOR_LMAX) then
         off = -1
      else if (kappa < 0) then
         off = OFF_GT(l)
      else
         off = OFF_LT(l)
      end if
   end subroutine coeff_offset

   ! ------------------------------------------------------------------
   ! The a_* kernels: split-real in, split-real out, as the drivers use
   ! them.  gx/gy/gz/g1 are the four Cartesian blocks a spin-included
   ! transform contracts; the spin-free ones read only g1.
   ! ------------------------------------------------------------------

   ! Bra side, spin included.  Writes the upper component at gsp(0:) and the
   ! lower at gsp(nket*nd:), which is how the C splits one buffer in two.
   subroutine a_bra_cart2spinor_si(gspR, gspI, gx, gy, gz, g1, nket, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: gx(0:), gy(0:), gz(0:), g1(0:)
      integer,  intent(in)    :: nket, kappa, l
      integer  :: nf, nd, ob, co, i, j, n
      real(dp) :: saR, saI, sbR, sbI, caR, caI, cbR, cbI, v1, vx, vy, vz
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ob = nket * nd
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
      do i = 0, nd - 1
         saR = 0.0_dp; saI = 0.0_dp; sbR = 0.0_dp; sbI = 0.0_dp
         do n = 0, nf - 1
            v1 = g1(j*nf+n); vx = gx(j*nf+n); vy = gy(j*nf+n); vz = gz(j*nf+n)
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            saR = saR + ( caR*v1 + caI*vz - cbR*vy + cbI*vx)
            saI = saI + (-caI*v1 + caR*vz + cbI*vy + cbR*vx)
            sbR = sbR + ( cbR*v1 - cbI*vz + caR*vy + caI*vx)
            sbI = sbI + (-cbI*v1 - cbR*vz - caI*vy + caR*vx)
         end do
         gspR(j*nd+i)      = saR
         gspI(j*nd+i)      = saI
         gspR(ob + j*nd+i) = sbR
         gspI(ob + j*nd+i) = sbI
      end do
      end do
   end subroutine a_bra_cart2spinor_si

   ! Bra side, spin free: only the scalar block contributes.
   subroutine a_bra_cart2spinor_sf(gspR, gspI, g1, nket, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: g1(0:)
      integer,  intent(in)    :: nket, kappa, l
      integer  :: nf, nd, ob, co, i, j, n
      real(dp) :: saR, saI, sbR, sbI, caR, caI, cbR, cbI, v1
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ob = nket * nd
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
      do i = 0, nd - 1
         saR = 0.0_dp; saI = 0.0_dp; sbR = 0.0_dp; sbI = 0.0_dp
         do n = 0, nf - 1
            v1 = g1(j*nf+n)
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            saR = saR + caR*v1
            saI = saI - caI*v1
            sbR = sbR + cbR*v1
            sbI = sbI - cbI*v1
         end do
         gspR(j*nd+i)      = saR
         gspI(j*nd+i)      = saI
         gspR(ob + j*nd+i) = sbR
         gspI(ob + j*nd+i) = sbI
      end do
      end do
   end subroutine a_bra_cart2spinor_sf

   ! Ket side, spin included.  lds is the leading dimension of the output,
   ! which is not nbra when the caller is placing a block into a larger array.
   subroutine a_ket_cart2spinor_si(gspR, gspI, gx, gy, gz, g1, lds, nbra, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: gx(0:), gy(0:), gz(0:), g1(0:)
      integer,  intent(in)    :: lds, nbra, kappa, l
      integer  :: nf, nd, ob, co, i, j, n
      real(dp) :: caR, caI, cbR, cbI, v1, vx, vy, vz
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ob = lds * nd
      call coeff_offset(kappa, l, co)

      do i = 0, nd - 1
         do j = 0, nbra - 1
            gspR(j+i*lds)      = 0.0_dp
            gspI(j+i*lds)      = 0.0_dp
            gspR(ob + j+i*lds) = 0.0_dp
            gspI(ob + j+i*lds) = 0.0_dp
         end do
         do n = 0, nf - 1
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            do j = 0, nbra - 1
               v1 = g1(j+n*nbra); vx = gx(j+n*nbra)
               vy = gy(j+n*nbra); vz = gz(j+n*nbra)
               ! [ 1+iz,  y+ix ] [ca]
               ! [-y+ix, 1-iz  ] [cb]
               gspR(j+i*lds)      = gspR(j+i*lds)      + ( caR*v1 - caI*vz + cbR*vy - cbI*vx)
               gspI(j+i*lds)      = gspI(j+i*lds)      + ( caI*v1 + caR*vz + cbI*vy + cbR*vx)
               gspR(ob + j+i*lds) = gspR(ob + j+i*lds) + ( cbR*v1 + cbI*vz - caR*vy - caI*vx)
               gspI(ob + j+i*lds) = gspI(ob + j+i*lds) + ( cbI*v1 - cbR*vz - caI*vy + caR*vx)
            end do
         end do
      end do
   end subroutine a_ket_cart2spinor_si

   subroutine a_ket_cart2spinor_sf(gspR, gspI, g1, lds, nbra, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: g1(0:)
      integer,  intent(in)    :: lds, nbra, kappa, l
      integer  :: nf, nd, ob, co, i, j, n
      real(dp) :: caR, caI, cbR, cbI, v1
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ob = lds * nd
      call coeff_offset(kappa, l, co)

      do i = 0, nd - 1
         do j = 0, nbra - 1
            gspR(j+i*lds)      = 0.0_dp
            gspI(j+i*lds)      = 0.0_dp
            gspR(ob + j+i*lds) = 0.0_dp
            gspI(ob + j+i*lds) = 0.0_dp
         end do
         do n = 0, nf - 1
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            do j = 0, nbra - 1
               v1 = g1(j+n*nbra)
               gspR(j+i*lds)      = gspR(j+i*lds)      + caR*v1
               gspI(j+i*lds)      = gspI(j+i*lds)      + caI*v1
               gspR(ob + j+i*lds) = gspR(ob + j+i*lds) + cbR*v1
               gspI(ob + j+i*lds) = gspI(ob + j+i*lds) + cbI*v1
            end do
         end do
      end do
   end subroutine a_ket_cart2spinor_sf

   ! ------------------------------------------------------------------
   ! The a_bra1_* family.  Same transforms as a_bra_*, but with a trailing
   ! index -- `ngrids` in the C, because int1e_grids introduced it -- that
   ! the 2e path reuses as the k and l block width.  The scalar versions
   ! (si, sf) read real Cartesian blocks; the z versions (zi, zf) read
   ! complex ones, held as two halves.
   ! ------------------------------------------------------------------

   subroutine a_bra1_cart2spinor_si(gspR, gspI, gx, gy, gz, g1, ngrids, nket, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: gx(0:), gy(0:), gz(0:), g1(0:)
      integer,  intent(in)    :: ngrids, nket, kappa, l
      integer  :: nf, nd, ndg, ob, co, i, j, n, m, s, t
      real(dp) :: caR, caI, cbR, cbI, v1, vx, vy, vz
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ndg = nd * ngrids
      ob = nket * ndg
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
         do i = 0, ndg - 1
            gspR(j*ndg+i)      = 0.0_dp
            gspI(j*ndg+i)      = 0.0_dp
            gspR(ob + j*ndg+i) = 0.0_dp
            gspI(ob + j*ndg+i) = 0.0_dp
         end do
         do i = 0, nd - 1
         do n = 0, nf - 1
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            s = (j*nf+n)*ngrids
            t = (j*nd+i)*ngrids
            do m = 0, ngrids - 1
               v1 = g1(s+m); vx = gx(s+m); vy = gy(s+m); vz = gz(s+m)
               gspR(t+m)      = gspR(t+m)      + ( caR*v1 + caI*vz - cbR*vy + cbI*vx)
               gspI(t+m)      = gspI(t+m)      + (-caI*v1 + caR*vz + cbI*vy + cbR*vx)
               gspR(ob + t+m) = gspR(ob + t+m) + ( cbR*v1 - cbI*vz + caR*vy + caI*vx)
               gspI(ob + t+m) = gspI(ob + t+m) + (-cbI*v1 - cbR*vz - caI*vy + caR*vx)
            end do
         end do
         end do
      end do
   end subroutine a_bra1_cart2spinor_si

   subroutine a_bra1_cart2spinor_sf(gspR, gspI, g1, ngrids, nket, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: g1(0:)
      integer,  intent(in)    :: ngrids, nket, kappa, l
      integer  :: nf, nd, ndg, ob, co, i, j, n, m, s, t
      real(dp) :: caR, caI, cbR, cbI, v1
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ndg = nd * ngrids
      ob = nket * ndg
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
         do i = 0, ndg - 1
            gspR(j*ndg+i)      = 0.0_dp
            gspI(j*ndg+i)      = 0.0_dp
            gspR(ob + j*ndg+i) = 0.0_dp
            gspI(ob + j*ndg+i) = 0.0_dp
         end do
         do i = 0, nd - 1
         do n = 0, nf - 1
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            s = (j*nf+n)*ngrids
            t = (j*nd+i)*ngrids
            do m = 0, ngrids - 1
               v1 = g1(s+m)
               gspR(t+m)      = gspR(t+m)      + caR*v1
               gspI(t+m)      = gspI(t+m)      - caI*v1
               gspR(ob + t+m) = gspR(ob + t+m) + cbR*v1
               gspI(ob + t+m) = gspI(ob + t+m) - cbI*v1
            end do
         end do
         end do
      end do
   end subroutine a_bra1_cart2spinor_sf

   ! Complex Cartesian in, spin included.  The 2x2 sigma matrix
   !   [ 1+iz,  y+ix ]
   !   [-y+ix, 1-iz  ]
   ! is applied first, then the conjugated coefficients.
   subroutine a_bra1_cart2spinor_zi(gspR, gspI, gx, gy, gz, g1, ngrids, nket, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: gx(0:), gy(0:), gz(0:), g1(0:)
      integer,  intent(in)    :: ngrids, nket, kappa, l
      integer  :: nf, nd, ndg, ob, oc, co, i, j, n, m, s, t
      real(dp) :: caR, caI, cbR, cbI
      real(dp) :: v11R, v12R, v21R, v22R, v11I, v12I, v21I, v22I
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ndg = nd * ngrids
      ob = nket * ndg
      oc = nket * nf * ngrids     ! start of the imaginary half of the input
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
         do i = 0, ndg - 1
            gspR(j*ndg+i)      = 0.0_dp
            gspI(j*ndg+i)      = 0.0_dp
            gspR(ob + j*ndg+i) = 0.0_dp
            gspI(ob + j*ndg+i) = 0.0_dp
         end do
         do i = 0, nd - 1
         do n = 0, nf - 1
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            s = (j*nf+n)*ngrids
            t = (j*nd+i)*ngrids
            do m = 0, ngrids - 1
               v11R = g1(s+m)      - gz(oc+s+m)
               v11I = g1(oc+s+m)   + gz(s+m)
               v12R = gy(s+m)      - gx(oc+s+m)
               v12I = gy(oc+s+m)   + gx(s+m)
               v21R = -gy(s+m)     - gx(oc+s+m)
               v21I = -gy(oc+s+m)  + gx(s+m)
               v22R = g1(s+m)      + gz(oc+s+m)
               v22I = g1(oc+s+m)   - gz(s+m)
               gspR(t+m)      = gspR(t+m)      + (caR*v11R + caI*v11I + cbR*v21R + cbI*v21I)
               gspI(t+m)      = gspI(t+m)      + (caR*v11I - caI*v11R + cbR*v21I - cbI*v21R)
               gspR(ob + t+m) = gspR(ob + t+m) + (caR*v12R + caI*v12I + cbR*v22R + cbI*v22I)
               gspI(ob + t+m) = gspI(ob + t+m) + (caR*v12I - caI*v12R + cbR*v22I - cbI*v22R)
            end do
         end do
         end do
      end do
   end subroutine a_bra1_cart2spinor_zi

   subroutine a_bra1_cart2spinor_zf(gspR, gspI, g1, ngrids, nket, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: g1(0:)
      integer,  intent(in)    :: ngrids, nket, kappa, l
      integer  :: nf, nd, ndg, ob, oc, co, i, j, n, m, s, t
      real(dp) :: caR, caI, cbR, cbI, v1R, v1I
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ndg = nd * ngrids
      ob = nket * ndg
      oc = nket * nf * ngrids
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
         do i = 0, ndg - 1
            gspR(j*ndg+i)      = 0.0_dp
            gspI(j*ndg+i)      = 0.0_dp
            gspR(ob + j*ndg+i) = 0.0_dp
            gspI(ob + j*ndg+i) = 0.0_dp
         end do
         do i = 0, nd - 1
         do n = 0, nf - 1
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            s = (j*nf+n)*ngrids
            t = (j*nd+i)*ngrids
            do m = 0, ngrids - 1
               v1R = g1(s+m); v1I = g1(oc+s+m)
               gspR(t+m)      = gspR(t+m)      + (caR*v1R + caI*v1I)
               gspI(t+m)      = gspI(t+m)      + (caR*v1I - caI*v1R)
               gspR(ob + t+m) = gspR(ob + t+m) + (cbR*v1R + cbI*v1I)
               gspI(ob + t+m) = gspI(ob + t+m) + (cbR*v1I - cbI*v1R)
            end do
         end do
         end do
      end do
   end subroutine a_bra1_cart2spinor_zf

   ! Ket side over an already-complex Cartesian block.  The C skips the
   ! multiply when a coefficient is zero; that is a speed test, not a
   ! numerical one -- adding 0*x changes nothing unless x is a NaN or an
   ! infinity, and neither reaches here -- so it is kept, in the same shape,
   ! rather than reasoned away.
   subroutine a_ket_cart2spinor(gspR, gspI, gcartR, gcartI, nbra, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: gcartR(0:), gcartI(0:)
      integer,  intent(in)    :: nbra, kappa, l
      integer  :: nf2, nd, co, i, j, n
      real(dp) :: cR, cI, gR, gI
      nf2 = len_cart(l) * 2
      nd = cint_len_spinor_kl(kappa, l)
      call coeff_offset(kappa, l, co)

      do i = 0, nd - 1
         do j = 0, nbra - 1
            gspR(j+i*nbra) = 0.0_dp
            gspI(j+i*nbra) = 0.0_dp
         end do
         do n = 0, nf2 - 1
            cR = g_trans_cart2jR(co + i*nf2 + n)
            cI = g_trans_cart2jI(co + i*nf2 + n)
            if (cR /= 0.0_dp) then
               if (cI /= 0.0_dp) then
                  do j = 0, nbra - 1
                     gR = gcartR(j+n*nbra); gI = gcartI(j+n*nbra)
                     gspR(j+i*nbra) = gspR(j+i*nbra) + (cR*gR - cI*gI)
                     gspI(j+i*nbra) = gspI(j+i*nbra) + (cI*gR + cR*gI)
                  end do
               else
                  do j = 0, nbra - 1
                     gR = gcartR(j+n*nbra); gI = gcartI(j+n*nbra)
                     gspR(j+i*nbra) = gspR(j+i*nbra) + cR*gR
                     gspI(j+i*nbra) = gspI(j+i*nbra) + cR*gI
                  end do
               end if
            else if (cI /= 0.0_dp) then
               do j = 0, nbra - 1
                  gR = gcartR(j+n*nbra); gI = gcartI(j+n*nbra)
                  gspR(j+i*nbra) = gspR(j+i*nbra) - cI*gI
                  gspI(j+i*nbra) = gspI(j+i*nbra) + cI*gR
               end do
            end if
         end do
      end do
   end subroutine a_ket_cart2spinor

   ! The same with a factor of i, which the C gets by swapping the two output
   ! halves and negating the real one.  Kept as that rather than as a
   ! multiplication, because it is exact and the multiplication is not
   ! obviously so once -0.0 is in play.
   subroutine a_iket_cart2spinor(gspR, gspI, gcartR, gcartI, nbra, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: gcartR(0:), gcartI(0:)
      integer,  intent(in)    :: nbra, kappa, l
      integer :: sz, i
      call a_ket_cart2spinor(gspI, gspR, gcartR, gcartI, nbra, kappa, l)
      sz = cint_len_spinor_kl(kappa, l) * nbra
      do i = 0, sz - 1
         gspR(i) = -gspR(i)
      end do
   end subroutine a_iket_cart2spinor

   subroutine a_ket1_cart2spinor(gspR, gspI, gcartR, gcartI, nbra, counts, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: gcartR(0:), gcartI(0:)
      integer,  intent(in)    :: nbra, counts, kappa, l
      integer  :: nf, nf2, nd, nds, nfs, ob, co, i, j, k, n
      real(dp) :: caR, caI, cbR, cbI, gaR, gaI, gbR, gbI
      nf = len_cart(l); nf2 = nf * 2
      nd = cint_len_spinor_kl(kappa, l)
      nds = nd * nbra
      nfs = nf * nbra
      ob = nfs * counts
      call coeff_offset(kappa, l, co)

      do i = 0, nd - 1
         do k = 0, counts - 1
         do j = 0, nbra - 1
            gspR(k*nds + j+i*nbra) = 0.0_dp
            gspI(k*nds + j+i*nbra) = 0.0_dp
         end do
         end do
         do n = 0, nf - 1
            caR = g_trans_cart2jR(co + i*nf2      + n)
            caI = g_trans_cart2jI(co + i*nf2      + n)
            cbR = g_trans_cart2jR(co + i*nf2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf2 + nf + n)
            do k = 0, counts - 1
            do j = 0, nbra - 1
               gaR = gcartR(k*nfs + j+n*nbra)
               gaI = gcartI(k*nfs + j+n*nbra)
               gbR = gcartR(ob + k*nfs + j+n*nbra)
               gbI = gcartI(ob + k*nfs + j+n*nbra)
               gspR(k*nds + j+i*nbra) = gspR(k*nds + j+i*nbra) &
                  + (caR*gaR - caI*gaI + cbR*gbR - cbI*gbI)
               gspI(k*nds + j+i*nbra) = gspI(k*nds + j+i*nbra) &
                  + (caR*gaI + caI*gaR + cbR*gbI + cbI*gbR)
            end do
            end do
         end do
      end do
   end subroutine a_ket1_cart2spinor

   subroutine a_iket1_cart2spinor(gspR, gspI, gcartR, gcartI, nbra, counts, kappa, l)
      real(dp), intent(inout) :: gspR(0:), gspI(0:)
      real(dp), intent(in)    :: gcartR(0:), gcartI(0:)
      integer,  intent(in)    :: nbra, counts, kappa, l
      integer :: sz, i
      call a_ket1_cart2spinor(gspI, gspR, gcartR, gcartI, nbra, counts, kappa, l)
      sz = cint_len_spinor_kl(kappa, l) * nbra * counts
      do i = 0, sz - 1
         gspR(i) = -gspR(i)
      end do
   end subroutine a_iket1_cart2spinor

   ! ------------------------------------------------------------------
   ! The CINTc2s_*_spinor* entry points.  These are what the generated
   ! relativistic integrals call, so they take and return complex(dp).
   ! ------------------------------------------------------------------

   ! Real Cartesian in, spinor out, spin free.  The bra transform conjugates:
   ! note the minus on every imaginary term.
   subroutine c2s_bra_spinor_e1sf(gsp, nket, gcart, kappa, l)
      complex(dp), intent(inout) :: gsp(0:)
      real(dp),    intent(in)    :: gcart(0:)
      integer,     intent(in)    :: nket, kappa, l
      integer  :: nf, nd, ob, co, i, j, n
      real(dp) :: saR, saI, sbR, sbI, caR, caI, cbR, cbI, v1
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ob = nket * nd
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
      do i = 0, nd - 1
         saR = 0.0_dp; saI = 0.0_dp; sbR = 0.0_dp; sbI = 0.0_dp
         do n = 0, nf - 1
            v1 = gcart(j*nf+n)
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            saR = saR + caR*v1
            saI = saI - caI*v1
            sbR = sbR + cbR*v1
            sbI = sbI - cbI*v1
         end do
         gsp(j*nd+i)      = cmplx(saR, saI, dp)
         gsp(ob + j*nd+i) = cmplx(sbR, sbI, dp)
      end do
      end do
   end subroutine c2s_bra_spinor_e1sf

   subroutine c2s_bra_spinor_sf(gsp, nket, gcart, kappa, l)
      complex(dp), intent(inout) :: gsp(0:)
      complex(dp), intent(in)    :: gcart(0:)
      integer,     intent(in)    :: nket, kappa, l
      integer  :: nf, nd, ob, co, i, j, n
      real(dp) :: saR, saI, sbR, sbI, caR, caI, cbR, cbI, v1R, v1I
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ob = nket * nd
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
      do i = 0, nd - 1
         saR = 0.0_dp; saI = 0.0_dp; sbR = 0.0_dp; sbI = 0.0_dp
         do n = 0, nf - 1
            v1R = real(gcart(j*nf+n), dp)
            v1I = aimag(gcart(j*nf+n))
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            saR = saR + (caR*v1R + caI*v1I)
            saI = saI + (caR*v1I - caI*v1R)
            sbR = sbR + (cbR*v1R + cbI*v1I)
            sbI = sbI + (cbR*v1I - cbI*v1R)
         end do
         gsp(j*nd+i)      = cmplx(saR, saI, dp)
         gsp(ob + j*nd+i) = cmplx(sbR, sbI, dp)
      end do
      end do
   end subroutine c2s_bra_spinor_sf

   subroutine c2s_bra_spinor(gsp, nket, gcart, kappa, l)
      complex(dp), intent(inout) :: gsp(0:)
      complex(dp), intent(in)    :: gcart(0:)
      integer,     intent(in)    :: nket, kappa, l
      integer  :: nf2, nd, co, i, j, n
      real(dp) :: sR, sI, cR, cI, gR, gI
      nf2 = len_cart(l) * 2
      nd = cint_len_spinor_kl(kappa, l)
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
      do i = 0, nd - 1
         sR = 0.0_dp; sI = 0.0_dp
         do n = 0, nf2 - 1
            gR = real(gcart(j*nf2+n), dp)
            gI = aimag(gcart(j*nf2+n))
            cR = g_trans_cart2jR(co + i*nf2 + n)
            cI = g_trans_cart2jI(co + i*nf2 + n)
            sR = sR + (cR*gR + cI*gI)
            sI = sI + (cR*gI - cI*gR)
         end do
         gsp(j*nd+i) = cmplx(sR, sI, dp)
      end do
      end do
   end subroutine c2s_bra_spinor

   ! Spin included: the Cartesian block is two halves, upper and lower.
   subroutine c2s_bra_spinor_si(gsp, nket, gcart, kappa, l)
      complex(dp), intent(inout) :: gsp(0:)
      complex(dp), intent(in)    :: gcart(0:)
      integer,     intent(in)    :: nket, kappa, l
      integer  :: nf, nd, ob, co, i, j, n
      real(dp) :: sR, sI, caR, caI, cbR, cbI, gaR, gaI, gbR, gbI
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ob = nf * nket
      call coeff_offset(kappa, l, co)

      do j = 0, nket - 1
      do i = 0, nd - 1
         sR = 0.0_dp; sI = 0.0_dp
         do n = 0, nf - 1
            caR = g_trans_cart2jR(co + i*nf*2      + n)
            caI = g_trans_cart2jI(co + i*nf*2      + n)
            cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
            cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
            gaR = real(gcart(j*nf+n), dp);      gaI = aimag(gcart(j*nf+n))
            gbR = real(gcart(ob + j*nf+n), dp); gbI = aimag(gcart(ob + j*nf+n))
            sR = sR + (caR*gaR + caI*gaI + cbR*gbR + cbI*gbI)
            sI = sI + (caR*gaI - caI*gaR + cbR*gbI - cbI*gbR)
         end do
         gsp(j*nd+i) = cmplx(sR, sI, dp)
      end do
      end do
   end subroutine c2s_bra_spinor_si

   subroutine c2s_ket_spinor(gsp, nbra, gcart, kappa, l)
      complex(dp), intent(inout) :: gsp(0:)
      complex(dp), intent(in)    :: gcart(0:)
      integer,     intent(in)    :: nbra, kappa, l
      integer  :: nf2, nd, co, i, j, n
      real(dp) :: cR, cI, gR, gI, aR, aI
      nf2 = len_cart(l) * 2
      nd = cint_len_spinor_kl(kappa, l)
      call coeff_offset(kappa, l, co)

      do i = 0, nd - 1
         do j = 0, nbra - 1
            gsp(j+i*nbra) = (0.0_dp, 0.0_dp)
         end do
         do n = 0, nf2 - 1
            cR = g_trans_cart2jR(co + i*nf2 + n)
            cI = g_trans_cart2jI(co + i*nf2 + n)
            do j = 0, nbra - 1
               gR = real(gcart(j+n*nbra), dp); gI = aimag(gcart(j+n*nbra))
               aR = real(gsp(j+i*nbra), dp);   aI = aimag(gsp(j+i*nbra))
               gsp(j+i*nbra) = cmplx(aR + (cR*gR - cI*gI), &
                                     aI + (cI*gR + cR*gI), dp)
            end do
         end do
      end do
   end subroutine c2s_ket_spinor

   subroutine c2s_iket_spinor(gsp, nbra, gcart, kappa, l)
      complex(dp), intent(inout) :: gsp(0:)
      complex(dp), intent(in)    :: gcart(0:)
      integer,     intent(in)    :: nbra, kappa, l
      integer  :: nf2, nd, co, i, j, n
      real(dp) :: cR, cI, gR, gI, aR, aI
      nf2 = len_cart(l) * 2
      nd = cint_len_spinor_kl(kappa, l)
      call coeff_offset(kappa, l, co)

      do i = 0, nd - 1
         do j = 0, nbra - 1
            gsp(j+i*nbra) = (0.0_dp, 0.0_dp)
         end do
         do n = 0, nf2 - 1
            cR = g_trans_cart2jR(co + i*nf2 + n)
            cI = g_trans_cart2jI(co + i*nf2 + n)
            do j = 0, nbra - 1
               gR = real(gcart(j+n*nbra), dp); gI = aimag(gcart(j+n*nbra))
               aR = real(gsp(j+i*nbra), dp);   aI = aimag(gsp(j+i*nbra))
               gsp(j+i*nbra) = cmplx(aR - (cI*gR + cR*gI), &
                                     aI + (cR*gR - cI*gI), dp)
            end do
         end do
      end do
   end subroutine c2s_iket_spinor

   ! gspa and gspb are the upper and lower components of the two-component
   ! vector, held in separate arrays rather than two halves of one.
   subroutine c2s_ket_spinor_sf1(gspa, gspb, gcart, lds, ldc, nctr, kappa, l)
      complex(dp), intent(inout) :: gspa(0:), gspb(0:)
      real(dp),    intent(in)    :: gcart(0:)
      integer,     intent(in)    :: lds, ldc, nctr, kappa, l
      integer  :: nf, nd, co, i, j, k, n, os, oc
      real(dp) :: caR, caI, cbR, cbI, v1, aR, aI, bR, bI
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      call coeff_offset(kappa, l, co)

      os = 0; oc = 0
      do k = 0, nctr - 1
         do i = 0, nd - 1
            do j = 0, ldc - 1
               gspa(os + j+i*lds) = (0.0_dp, 0.0_dp)
               gspb(os + j+i*lds) = (0.0_dp, 0.0_dp)
            end do
            do n = 0, nf - 1
               caR = g_trans_cart2jR(co + i*nf*2      + n)
               caI = g_trans_cart2jI(co + i*nf*2      + n)
               cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
               cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
               do j = 0, ldc - 1
                  v1 = gcart(oc + j+n*ldc)
                  aR = real(gspa(os + j+i*lds), dp); aI = aimag(gspa(os + j+i*lds))
                  bR = real(gspb(os + j+i*lds), dp); bI = aimag(gspb(os + j+i*lds))
                  gspa(os + j+i*lds) = cmplx(aR + caR*v1, aI + caI*v1, dp)
                  gspb(os + j+i*lds) = cmplx(bR + cbR*v1, bI + cbI*v1, dp)
               end do
            end do
         end do
         os = os + nd * lds
         oc = oc + nf * ldc
      end do
   end subroutine c2s_ket_spinor_sf1

   subroutine c2s_iket_spinor_sf1(gspa, gspb, gcart, lds, ldc, nctr, kappa, l)
      complex(dp), intent(inout) :: gspa(0:), gspb(0:)
      real(dp),    intent(in)    :: gcart(0:)
      integer,     intent(in)    :: lds, ldc, nctr, kappa, l
      integer  :: nf, nd, co, i, j, k, n, os, oc
      real(dp) :: caR, caI, cbR, cbI, v1, aR, aI, bR, bI
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      call coeff_offset(kappa, l, co)

      os = 0; oc = 0
      do k = 0, nctr - 1
         do i = 0, nd - 1
            do j = 0, ldc - 1
               gspa(os + j+i*lds) = (0.0_dp, 0.0_dp)
               gspb(os + j+i*lds) = (0.0_dp, 0.0_dp)
            end do
            do n = 0, nf - 1
               caR = g_trans_cart2jR(co + i*nf*2      + n)
               caI = g_trans_cart2jI(co + i*nf*2      + n)
               cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
               cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
               do j = 0, ldc - 1
                  v1 = gcart(oc + j+n*ldc)
                  aR = real(gspa(os + j+i*lds), dp); aI = aimag(gspa(os + j+i*lds))
                  bR = real(gspb(os + j+i*lds), dp); bI = aimag(gspb(os + j+i*lds))
                  gspa(os + j+i*lds) = cmplx(aR - caI*v1, aI + caR*v1, dp)
                  gspb(os + j+i*lds) = cmplx(bR - cbI*v1, bI + cbR*v1, dp)
               end do
            end do
         end do
         os = os + nd * lds
         oc = oc + nf * ldc
      end do
   end subroutine c2s_iket_spinor_sf1

   ! Spin included: gcart holds four blocks in the order x, y, z, 1 -- the
   ! sigma components and then the scalar.
   subroutine c2s_ket_spinor_si1(gspa, gspb, gcart, lds, ldc, nctr, kappa, l)
      complex(dp), intent(inout) :: gspa(0:), gspb(0:)
      real(dp),    intent(in)    :: gcart(0:)
      integer,     intent(in)    :: lds, ldc, nctr, kappa, l
      integer  :: nf, nd, ngc, co, i, j, k, n, os, ox, oy, oz, o1
      real(dp) :: caR, caI, cbR, cbI, v1, vx, vy, vz, aR, aI, bR, bI
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ngc = nf * ldc
      call coeff_offset(kappa, l, co)

      os = 0
      ox = 0; oy = nctr*ngc; oz = 2*nctr*ngc; o1 = 3*nctr*ngc
      do k = 0, nctr - 1
         do i = 0, nd - 1
            do j = 0, ldc - 1
               gspa(os + j+i*lds) = (0.0_dp, 0.0_dp)
               gspb(os + j+i*lds) = (0.0_dp, 0.0_dp)
            end do
            do n = 0, nf - 1
               caR = g_trans_cart2jR(co + i*nf*2      + n)
               caI = g_trans_cart2jI(co + i*nf*2      + n)
               cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
               cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
               do j = 0, ldc - 1
                  v1 = gcart(o1 + j+n*ldc); vx = gcart(ox + j+n*ldc)
                  vy = gcart(oy + j+n*ldc); vz = gcart(oz + j+n*ldc)
                  aR = real(gspa(os + j+i*lds), dp); aI = aimag(gspa(os + j+i*lds))
                  bR = real(gspb(os + j+i*lds), dp); bI = aimag(gspb(os + j+i*lds))
                  gspa(os + j+i*lds) = cmplx( &
                     aR + ( caR*v1 - caI*vz + cbR*vy - cbI*vx), &
                     aI + ( caI*v1 + caR*vz + cbI*vy + cbR*vx), dp)
                  gspb(os + j+i*lds) = cmplx( &
                     bR + ( cbR*v1 + cbI*vz - caR*vy - caI*vx), &
                     bI + ( cbI*v1 - cbR*vz - caI*vy + caR*vx), dp)
               end do
            end do
         end do
         os = os + nd * lds
         ox = ox + ngc; oy = oy + ngc; oz = oz + ngc; o1 = o1 + ngc
      end do
   end subroutine c2s_ket_spinor_si1

   subroutine c2s_iket_spinor_si1(gspa, gspb, gcart, lds, ldc, nctr, kappa, l)
      complex(dp), intent(inout) :: gspa(0:), gspb(0:)
      real(dp),    intent(in)    :: gcart(0:)
      integer,     intent(in)    :: lds, ldc, nctr, kappa, l
      integer  :: nf, nd, ngc, co, i, j, k, n, os, ox, oy, oz, o1
      real(dp) :: caR, caI, cbR, cbI, v1, vx, vy, vz, aR, aI, bR, bI
      nf = len_cart(l)
      nd = cint_len_spinor_kl(kappa, l)
      ngc = nf * ldc
      call coeff_offset(kappa, l, co)

      os = 0
      ox = 0; oy = nctr*ngc; oz = 2*nctr*ngc; o1 = 3*nctr*ngc
      do k = 0, nctr - 1
         do i = 0, nd - 1
            do j = 0, ldc - 1
               gspa(os + j+i*lds) = (0.0_dp, 0.0_dp)
               gspb(os + j+i*lds) = (0.0_dp, 0.0_dp)
            end do
            do n = 0, nf - 1
               caR = g_trans_cart2jR(co + i*nf*2      + n)
               caI = g_trans_cart2jI(co + i*nf*2      + n)
               cbR = g_trans_cart2jR(co + i*nf*2 + nf + n)
               cbI = g_trans_cart2jI(co + i*nf*2 + nf + n)
               do j = 0, ldc - 1
                  v1 = gcart(o1 + j+n*ldc); vx = gcart(ox + j+n*ldc)
                  vy = gcart(oy + j+n*ldc); vz = gcart(oz + j+n*ldc)
                  aR = real(gspa(os + j+i*lds), dp); aI = aimag(gspa(os + j+i*lds))
                  bR = real(gspb(os + j+i*lds), dp); bI = aimag(gspb(os + j+i*lds))
                  gspa(os + j+i*lds) = cmplx( &
                     aR - ( caI*v1 + caR*vz + cbI*vy + cbR*vx), &
                     aI + ( caR*v1 - caI*vz + cbR*vy - cbI*vx), dp)
                  gspb(os + j+i*lds) = cmplx( &
                     bR - ( cbI*v1 - cbR*vz - caI*vy + caR*vx), &
                     bI + ( cbR*v1 + cbI*vz - caR*vy - caI*vx), dp)
               end do
            end do
         end do
         os = os + nd * lds
         ox = ox + ngc; oy = oy + ngc; oz = oz + ngc; o1 = o1 + ngc
      end do
   end subroutine c2s_iket_spinor_si1

   ! ------------------------------------------------------------------

   ! Reorder the (i,k,l,j) block a 2e integral produces into the caller's
   ! (i,j,k,l) layout, merging the split real and imaginary halves as it goes.
   pure subroutine zcopy_iklj(fijkl, foff, gR, gI, ni, nj, nk, mi, mj, mk, ml)
      complex(dp), intent(inout) :: fijkl(0:)
      real(dp),    intent(in)    :: gR(0:), gI(0:)
      integer,     intent(in)    :: foff, ni, nj, nk, mi, mj, mk, ml
      integer :: i, j, k, l, nij, nijk, mik, mikl, fb, gb
      nij = ni*nj; nijk = nij*nk; mik = mi*mk; mikl = mik*ml
      fb = foff; gb = 0
      do l = 0, ml - 1
         do k = 0, mk - 1
            do j = 0, mj - 1
               do i = 0, mi - 1
                  fijkl(fb + k*nij + ni*j + i) = &
                     cmplx(gR(gb + k*mi + mikl*j + i), &
                           gI(gb + k*mi + mikl*j + i), dp)
               end do
            end do
         end do
         fb = fb + nijk
         gb = gb + mik
      end do
   end subroutine zcopy_iklj

   pure subroutine c2s_zset0(out, dims, counts)
      complex(dp), intent(inout) :: out(0:)
      integer,     intent(in)    :: dims(0:), counts(0:)
      integer :: i, j, k, l, ni, nij, nijk
      ni = dims(0); nij = ni*dims(1); nijk = nij*dims(2)
      do l = 0, counts(3) - 1
      do k = 0, counts(2) - 1
      do j = 0, counts(1) - 1
      do i = 0, counts(0) - 1
         out(l*nijk + k*nij + j*ni + i) = (0.0_dp, 0.0_dp)
      end do
      end do
      end do
      end do
   end subroutine c2s_zset0

end module cint_cart2spinor
