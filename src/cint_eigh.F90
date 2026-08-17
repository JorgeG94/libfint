!
! Symmetric tridiagonal eigenproblem for the Wheeler quadrature.
!
! This replaces src/eigh.c, and it is now a thin choice between two ways of
! doing that.  By default it calls cint_eigh_dstemr, which is eigh.c
! translated -- so the port needs no LAPACK at all, exactly as the C needs
! none, and is bit-identical here as everywhere else.
!
! WITH_EXTERNAL_LAPACK swaps in a real dstemr instead.  That is a supported
! configuration rather than a fallback: a caller that already links LAPACK
! (metalquicha does, through pic_lapack_interfaces) may prefer a vendor
! implementation, which will be faster and may be more accurate.  It is not
! the default because it is not bit-identical -- same algorithm, different
! implementation -- and the catalogue check then needs a tolerance above 5
! Rys roots.
!
module cint_eigh
   use cint_const, only: dp
   use cint_eigh_dstemr, only: dstemr_diagonalize
   implicit none
   private

   public :: cint_diagonalize

   integer, parameter :: MXRYSROOTS = 32

#ifdef WITH_EXTERNAL_LAPACK
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
#endif

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

#ifdef WITH_EXTERNAL_LAPACK
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
#else
      info = 0
      if (n <= 0) return
      ! the C reads n entries of the off-diagonal though only n-1 mean
      ! anything, and its own eigh.c zeroes the last one on the way in
      offd(n-1) = 0.0_dp
      info = dstemr_diagonalize(n, diag, offd, eig, vec)
#endif
   end function cint_diagonalize

end module cint_eigh
