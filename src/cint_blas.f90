!
! Small dense-array helpers.
!
! Replaces src/fblas.c, which despite the name calls no BLAS: it hand-rolls
! four small matrix multiplies and two transposes, because the blocks involved
! are shell-sized.
!
! WHY NOT pic_blas (or any BLAS) HERE.  Measured on this machine at the sizes
! that actually occur -- (m,n,k) of (6,5,6), (10,7,10), (15,9,15), (28,13,28),
! being l = 2, 3, 4 and 6 -- a plain Fortran loop against reference DGEMM:
!
!     flags                              6x5x6   10x7x10  15x9x15  28x13x28
!     -O2                                 0.95     0.66     0.47     0.37
!     -O3                                 1.25     1.05     0.91     0.78
!     -O3 -march=native -funroll-loops     1.31     1.25     1.08     1.04
!
! (ratio dgemm/loop; above 1 the loop wins).  So it is close, and which side
! wins depends on the flags.  More to the point, the matrix-multiply path is
! only reached for l >= 5: libcint hand-unrolls l = 0..4, which is every basis
! set metalquicha uses.  A BLAS dependency for a path that never executes is a
! poor trade, so the loops stay -- they live in cint_cart2sph, where they are
! used.  If that ever changes, pic_gemm is a drop-in at the two call sites.
!
! The transposes below are permutations, so the C's four-way unrolling is a
! scheduling detail with no effect on the values; the simple loops here are
! bit-identical to it.
!
module cint_blas
   use cint_const, only: dp
   implicit none
   private

   public :: cint_dset0, cint_daxpy2v, cint_dmat_transpose, cint_dplus_transpose

contains

   pure subroutine cint_dset0(n, x)
      integer,  intent(in)  :: n
      real(dp), intent(out) :: x(0:)
      integer :: i
      do i = 0, n - 1
         x(i) = 0.0_dp
      end do
   end subroutine cint_dset0

   ! v = a*x + y
   pure subroutine cint_daxpy2v(n, a, x, y, v)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: a, x(0:), y(0:)
      real(dp), intent(out) :: v(0:)
      integer :: i
      do i = 0, n - 1
         v(i) = a * x(i) + y(i)
      end do
   end subroutine cint_daxpy2v

   ! a(m,n) -> a_t(n,m), both flat and row-major as the C has them.
   pure subroutine cint_dmat_transpose(a_t, a, m, n)
      real(dp), intent(out) :: a_t(0:)
      real(dp), intent(in)  :: a(0:)
      integer,  intent(in)  :: m, n
      integer :: i, j
      do j = 0, n - 1
         do i = 0, m - 1
            a_t(j*m + i) = a(i*n + j)
         end do
      end do
   end subroutine cint_dmat_transpose

   ! As above but accumulating, which is how the multi-component drivers fold
   ! a tensor component into the running total.
   pure subroutine cint_dplus_transpose(a_t, a, m, n)
      real(dp), intent(inout) :: a_t(0:)
      real(dp), intent(in)    :: a(0:)
      integer,  intent(in)    :: m, n
      integer :: i, j
      do j = 0, n - 1
         do i = 0, m - 1
            a_t(j*m + i) = a_t(j*m + i) + a(i*n + j)
         end do
      end do
   end subroutine cint_dplus_transpose

end module cint_blas
