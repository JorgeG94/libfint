! Guarded by COMPILE_F12.  Without it this file compiles to nothing at all,
! which is what lets every source in fortran/src be listed unconditionally --
! and the reason it has to be a guard rather than "just always build it" is
! that the chain below ends at cint_tab_stg_roots, a 100 MB generated table
! that is deliberately not in the repository.  See
! scripts/stg_roots_to_fortran.py.
!
! No #else stub: fpm scans for module names lexically, before the
! preprocessor runs, so a second `module cint_stg_quad` in this file reads to it
! as a duplicate definition.  An absent module is also the better error for
! anyone who USEs it without asking for F12.
#ifdef COMPILE_F12
!
! The Slater-type-geminal quadrature.
!
! Ported from src/stg_roots.c, which libcint in turn took from libslater
! (https://github.com/nubakery/libslater).  Only built with WITH_F12, as in
! the C.
!
! The roots come out of a two-dimensional Chebyshev interpolation -- in
! log(t) and log10(u) -- over a table of 3,567,200 coefficients.  That table
! is generated rather than committed; scripts/stg_roots_to_fortran.py says
! why, and what it costs.
!
! The C unrolls all three kernels by hand: _clenshaw_dc across fourteen
! lanes in groups of four, _clenshaw_d1 across two roots at a time, and
! _matmul_14_14 across all fourteen outputs.  None of that unrolling changes
! the arithmetic within a lane, and the lanes are independent, so the loops
! here are bit-identical to it.  Written as loops because fourteen copies of
! the same line is how a transcription error hides.
!
module cint_stg_quad
   use cint_const, only: dp
   use cint_envs,  only: MXRYSROOTS
   use cint_tab_stg_roots, only: STG_X, STG_W
   implicit none
   private

   public :: cint_stg_roots

   ! cos(pi*j*(2*i+1)/28) for i, j in 0..13, the C's COS_14_14.
   real(dp), parameter :: COS_14_14(0:195) = [ &
      1.0_dp, 9.9371220989324260398e-01_dp, 9.7492791218182361934e-01_dp, &
      9.4388333030836757409e-01_dp, 9.0096886790241914600e-01_dp, 8.4672419922828412453e-01_dp, &
      7.8183148246802980363e-01_dp, 7.0710678118654757274e-01_dp, 6.2348980185873348336e-01_dp, &
      5.3203207651533657163e-01_dp, 4.3388373911755812040e-01_dp, 3.3027906195516709698e-01_dp, &
      2.2252093395631439288e-01_dp, 1.1196447610330785560e-01_dp, 1.0_dp, &
      9.4388333030836757409e-01_dp, 7.8183148246802980363e-01_dp, 5.3203207651533657163e-01_dp, &
      2.2252093395631439288e-01_dp, -1.1196447610330785560e-01_dp, -4.3388373911755812040e-01_dp, &
      -7.0710678118654757274e-01_dp, -9.0096886790241914600e-01_dp, -9.9371220989324260398e-01_dp, &
      -9.7492791218182361934e-01_dp, -8.4672419922828412453e-01_dp, -6.2348980185873348336e-01_dp, &
      -3.3027906195516709698e-01_dp, 1.0_dp, 8.4672419922828412453e-01_dp, &
      4.3388373911755812040e-01_dp, -1.1196447610330785560e-01_dp, -6.2348980185873348336e-01_dp, &
      -9.4388333030836757409e-01_dp, -9.7492791218182361934e-01_dp, -7.0710678118654757274e-01_dp, &
      -2.2252093395631439288e-01_dp, 3.3027906195516709698e-01_dp, 7.8183148246802980363e-01_dp, &
      9.9371220989324260398e-01_dp, 9.0096886790241914600e-01_dp, 5.3203207651533657163e-01_dp, &
      1.0_dp, 7.0710678118654757274e-01_dp, 0.0_dp, &
      -7.0710678118654757274e-01_dp, -1.0_dp, -7.0710678118654757274e-01_dp, &
      0.0_dp, 7.0710678118654757274e-01_dp, 1.0_dp, &
      7.0710678118654757274e-01_dp, 0.0_dp, -7.0710678118654757274e-01_dp, &
      -1.0_dp, -7.0710678118654757274e-01_dp, 1.0_dp, &
      5.3203207651533657163e-01_dp, -4.3388373911755812040e-01_dp, -9.9371220989324260398e-01_dp, &
      -6.2348980185873348336e-01_dp, 3.3027906195516709698e-01_dp, 9.7492791218182361934e-01_dp, &
      7.0710678118654757274e-01_dp, -2.2252093395631439288e-01_dp, -9.4388333030836757409e-01_dp, &
      -7.8183148246802980363e-01_dp, 1.1196447610330785560e-01_dp, 9.0096886790241914600e-01_dp, &
      8.4672419922828412453e-01_dp, 1.0_dp, 3.3027906195516709698e-01_dp, &
      -7.8183148246802980363e-01_dp, -8.4672419922828412453e-01_dp, 2.2252093395631439288e-01_dp, &
      9.9371220989324260398e-01_dp, 4.3388373911755812040e-01_dp, -7.0710678118654757274e-01_dp, &
      -9.0096886790241914600e-01_dp, 1.1196447610330785560e-01_dp, 9.7492791218182361934e-01_dp, &
      5.3203207651533657163e-01_dp, -6.2348980185873348336e-01_dp, -9.4388333030836757409e-01_dp, &
      1.0_dp, 1.1196447610330785560e-01_dp, -9.7492791218182361934e-01_dp, &
      -3.3027906195516709698e-01_dp, 9.0096886790241914600e-01_dp, 5.3203207651533657163e-01_dp, &
      -7.8183148246802980363e-01_dp, -7.0710678118654757274e-01_dp, 6.2348980185873348336e-01_dp, &
      8.4672419922828412453e-01_dp, -4.3388373911755812040e-01_dp, -9.4388333030836757409e-01_dp, &
      2.2252093395631439288e-01_dp, 9.9371220989324260398e-01_dp, 1.0_dp, &
      -1.1196447610330785560e-01_dp, -9.7492791218182361934e-01_dp, 3.3027906195516709698e-01_dp, &
      9.0096886790241914600e-01_dp, -5.3203207651533657163e-01_dp, -7.8183148246802980363e-01_dp, &
      7.0710678118654757274e-01_dp, 6.2348980185873348336e-01_dp, -8.4672419922828412453e-01_dp, &
      -4.3388373911755812040e-01_dp, 9.4388333030836757409e-01_dp, 2.2252093395631439288e-01_dp, &
      -9.9371220989324260398e-01_dp, 1.0_dp, -3.3027906195516709698e-01_dp, &
      -7.8183148246802980363e-01_dp, 8.4672419922828412453e-01_dp, 2.2252093395631439288e-01_dp, &
      -9.9371220989324260398e-01_dp, 4.3388373911755812040e-01_dp, 7.0710678118654757274e-01_dp, &
      -9.0096886790241914600e-01_dp, -1.1196447610330785560e-01_dp, 9.7492791218182361934e-01_dp, &
      -5.3203207651533657163e-01_dp, -6.2348980185873348336e-01_dp, 9.4388333030836757409e-01_dp, &
      1.0_dp, -5.3203207651533657163e-01_dp, -4.3388373911755812040e-01_dp, &
      9.9371220989324260398e-01_dp, -6.2348980185873348336e-01_dp, -3.3027906195516709698e-01_dp, &
      9.7492791218182361934e-01_dp, -7.0710678118654757274e-01_dp, -2.2252093395631439288e-01_dp, &
      9.4388333030836757409e-01_dp, -7.8183148246802980363e-01_dp, -1.1196447610330785560e-01_dp, &
      9.0096886790241914600e-01_dp, -8.4672419922828412453e-01_dp, 1.0_dp, &
      -7.0710678118654757274e-01_dp, 0.0_dp, 7.0710678118654757274e-01_dp, &
      -1.0_dp, 7.0710678118654757274e-01_dp, 0.0_dp, &
      -7.0710678118654757274e-01_dp, 1.0_dp, -7.0710678118654757274e-01_dp, &
      0.0_dp, 7.0710678118654757274e-01_dp, -1.0_dp, &
      7.0710678118654757274e-01_dp, 1.0_dp, -8.4672419922828412453e-01_dp, &
      4.3388373911755812040e-01_dp, 1.1196447610330785560e-01_dp, -6.2348980185873348336e-01_dp, &
      9.4388333030836757409e-01_dp, -9.7492791218182361934e-01_dp, 7.0710678118654757274e-01_dp, &
      -2.2252093395631439288e-01_dp, -3.3027906195516709698e-01_dp, 7.8183148246802980363e-01_dp, &
      -9.9371220989324260398e-01_dp, 9.0096886790241914600e-01_dp, -5.3203207651533657163e-01_dp, &
      1.0_dp, -9.4388333030836757409e-01_dp, 7.8183148246802980363e-01_dp, &
      -5.3203207651533657163e-01_dp, 2.2252093395631439288e-01_dp, 1.1196447610330785560e-01_dp, &
      -4.3388373911755812040e-01_dp, 7.0710678118654757274e-01_dp, -9.0096886790241914600e-01_dp, &
      9.9371220989324260398e-01_dp, -9.7492791218182361934e-01_dp, 8.4672419922828412453e-01_dp, &
      -6.2348980185873348336e-01_dp, 3.3027906195516709698e-01_dp, 1.0_dp, &
      -9.9371220989324260398e-01_dp, 9.7492791218182361934e-01_dp, -9.4388333030836757409e-01_dp, &
      9.0096886790241914600e-01_dp, -8.4672419922828412453e-01_dp, 7.8183148246802980363e-01_dp, &
      -7.0710678118654757274e-01_dp, 6.2348980185873348336e-01_dp, -5.3203207651533657163e-01_dp, &
      4.3388373911755812040e-01_dp, -3.3027906195516709698e-01_dp, 2.2252093395631439288e-01_dp, &
      -1.1196447610330785560e-01_dp ]

contains

   ! Clenshaw in the u direction: fourteen independent Chebyshev sums per
   ! root, each over fourteen coefficients with stride 14.
   pure subroutine clenshaw_dc(rr, x, xoff, u, nroot)
      real(dp), intent(out) :: rr(0:)
      real(dp), intent(in)  :: x(0:)
      integer,  intent(in)  :: xoff, nroot
      real(dp), intent(in)  :: u
      real(dp) :: d, g, u2
      integer  :: i, j, k, b
      u2 = u * 2.0_dp
      do i = 0, nroot - 1
         b = xoff + 196*i
         do j = 0, 13
            d = 0.0_dp
            g = x(b + 13 + 14*j)
            do k = 11, 1, -2
               d = u2 * g - d + x(b + k + 1 + 14*j)
               g = u2 * d - g + x(b + k     + 14*j)
            end do
            rr(j + 14*i) = u * g - d + x(b + 14*j) * 0.5_dp
         end do
      end do
   end subroutine clenshaw_dc

   ! The 14x14 cosine transform between the two Chebyshev directions.
   pure subroutine matmul_14_14(imc, im, nroot)
      real(dp), intent(out) :: imc(0:)
      real(dp), intent(in)  :: im(0:)
      integer,  intent(in)  :: nroot
      real(dp), parameter :: O7 = 0.14285714285714285714_dp
      real(dp) :: d0(0:13), s
      integer  :: i, j, m
      do i = 0, nroot - 1
         d0 = 0.0_dp
         do j = 0, 13
            s = im(j + 14*i)
            do m = 0, 13
               d0(m) = d0(m) + s * COS_14_14(j*14 + m)
            end do
         end do
         do m = 0, 13
            imc(m + 14*i) = O7 * d0(m)
         end do
      end do
   end subroutine matmul_14_14

   ! Clenshaw in the t direction: one Chebyshev sum per root.
   pure subroutine clenshaw_d1(rr, x, u, nroot)
      real(dp), intent(out) :: rr(0:)
      real(dp), intent(in)  :: x(0:)
      integer,  intent(in)  :: nroot
      real(dp), intent(in)  :: u
      real(dp) :: d0, g0, u2
      integer  :: i, k
      u2 = u * 2.0_dp
      do i = 0, nroot - 1
         d0 = 0.0_dp
         g0 = x(13 + 14*i)
         do k = 12, 1, -2
            d0 = u2 * g0 - d0 + x(k     + 14*i)
            g0 = u2 * d0 - g0 + x(k - 1 + 14*i)
         end do
         rr(i) = u * g0 - d0 + x(14*i) * 0.5_dp
      end do
   end subroutine clenshaw_d1

   subroutine cint_stg_roots(nroots, ta, ua, rr, ww)
      integer,  intent(in)    :: nroots
      real(dp), intent(in)    :: ta, ua
      real(dp), intent(inout) :: rr(0:), ww(0:)

      real(dp) :: u, uu, t, tt
      integer  :: k, iu, it, base, offset
      ! Fixed bound, not 14*nroots: an automatic array with a runtime bound
      ! goes on the heap under gfortran, and a malloc/free pair per call is
      ! what this costs.  nroots never exceeds MXRYSROOTS.
      real(dp) :: im(0:14*MXRYSROOTS-1), imc(0:14*MXRYSROOTS-1)

      base = (nroots-1)*nroots/2 * 19600

      t = ta
      if (t > 19682.99_dp) t = 19682.99_dp
      u = ua
      if (t > 1.0_dp) then
         tt = log(t) * 0.9102392266268373_dp + 1.0_dp    ! log(t)/log(3) + 1
      else
         tt = sqrt(t)
      end if
      uu = log10(u)

      it = int(tt)
      tt = tt - it
      tt = 2.0_dp * tt - 1.0_dp

      iu = int(uu + 7.0_dp)
      if (iu < 0 .or. iu > 10) then
         ! The C prints and exits here.  Same refusal, said the Fortran way.
         error stop "cint_stg_roots: the tabulation assumes 1.0e-7 < U < 1.0e3"
      end if
      uu = uu - (iu - 7)
      uu = 2.0_dp * uu - 1.0_dp

      offset = nroots * 196 * (iu + it * 10)
      call clenshaw_dc(im,  STG_X, base + offset, uu, nroots)
      call matmul_14_14(imc, im, nroots)
      call clenshaw_d1(rr, imc, tt, nroots)
      call clenshaw_dc(im,  STG_W, base + offset, uu, nroots)
      call matmul_14_14(imc, im, nroots)
      call clenshaw_d1(ww, imc, tt, nroots)

      uu = 1.0_dp / sqrt(ua)
      do k = 0, nroots - 1
         ww(k) = ww(k) * uu
      end do
   end subroutine cint_stg_roots

end module cint_stg_quad
#endif
