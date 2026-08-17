!
! Symmetric tridiagonal eigenproblem for the Wheeler quadrature.
!
! This replaces src/eigh.c, and replacing it is the point.  That file is 1,477
! lines, and reading its header explains why: it is LAPACK's DSTEMR, vendored
! and hand-translated, used only when the build cannot find a real LAPACK.
! When LAPACK_FOUND is set the C calls dstemr_ and the 1,477 lines are dead.
! libcint's CMakeLists never sets it, so the shipped library always runs the
! copy.
!
! metalquicha already links LAPACK through pic_lapack_interfaces, so calling
! dstemr directly costs nothing and deletes the largest single file in the
! port's dependency closure.
!
! ⚠️ This is the one place the port is deliberately not bit-identical to the C.
! Same algorithm, but a different implementation of it, so agreement is to
! within the eigensolver's own accuracy rather than to the last bit.  The
! Wheeler tests therefore check invariants -- that the rule integrates its own
! moments -- as well as comparing against the C.
!
module cint_eigh
   use cint_const, only: dp
   implicit none
   private

   public :: cint_diagonalize

   integer, parameter :: MXRYSROOTS = 32

   interface
      subroutine dstemr(jobz, range, n, d, e, vl, vu, il, iu, m, w, z, ldz, &
                        nzc, isuppz, tryrac, work, lwork, iwork, liwork, info)
         import :: dp
         character,  intent(in)    :: jobz, range
         integer,    intent(in)    :: n, il, iu, ldz, nzc, lwork, liwork
         real(dp),   intent(inout) :: d(*), e(*)
         real(dp),   intent(in)    :: vl, vu
         integer,    intent(out)   :: m, isuppz(*), iwork(*), info
         real(dp),   intent(out)   :: w(*), z(ldz, *), work(*)
         logical,    intent(inout) :: tryrac
      end subroutine dstemr
   end interface

contains

   ! Eigenvalues and eigenvectors of the tridiagonal matrix with diagonal
   ! diag(0:n-1) and off-diagonal offd(0:n-2).  Both are overwritten, as in the
   ! C.  `vec` comes back flat and 0-based with column-major LDZ = n, so the
   ! first component of eigenvector i is vec(i*n) -- which is what the weight
   ! formula wants.
   function cint_diagonalize(n, diag, offd, eig, vec) result(info)
      integer,  intent(in)    :: n
      real(dp), intent(inout) :: diag(0:), offd(0:)
      real(dp), intent(out)   :: eig(0:), vec(0:)
      integer :: info

      real(dp) :: vl, vu
      integer  :: il, iu, m
      integer  :: isuppz(2*MXRYSROOTS)
      logical  :: tryrac
      real(dp) :: work(18*MXRYSROOTS)
      integer  :: iwork(10*MXRYSROOTS)

      vl = 0.0_dp
      vu = 0.0_dp
      il = 0
      iu = 0
      m = 0
      tryrac = .true.
      info = 0

      if (n <= 0) return

      ! DSTEMR reads n entries of e even though only n-1 are meaningful.
      offd(n-1) = 0.0_dp

      call dstemr('V', 'A', n, diag, offd, vl, vu, il, iu, m, &
                  eig, vec, n, n, isuppz, tryrac, &
                  work, size(work), iwork, size(iwork), info)
   end function cint_diagonalize

end module cint_eigh
