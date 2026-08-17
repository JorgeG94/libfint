!
! Real roots of the Rys characteristic polynomial.
!
! Ported from src/find_roots.c, which is itself a cut-down mpmath
! matrices.eigen.eig restricted to real eigenvalues.  Two strategies: build the
! companion matrix and run an unshifted-deflation Hessenberg QR on it, and if
! that fails, fall back to bracketed root polishing on successive orders.
!
! Double precision only, deliberately.  The C's long double and __float128
! Schmidt paths both call this with a plain `double *` array, so there is no
! extended-precision variant to mirror.
!
! Indexing follows the C exactly: every array is 0-based and the 2D ones stay
! flat with the same row-major arithmetic (PORT_TO_FORTRAN.md 3.6).  `cs` is
! conceptually (nroots+1) x (nroots+1) and addressed as cs(order*nroots1 + i).
!
module cint_find_roots
   use cint_const, only: dp
   use iso_fortran_env, only: error_unit
   implicit none
   private

   public :: cint_polynomial_roots, poly_value1

   integer, parameter :: MXRYSROOTS = 32

contains

   ! Horner evaluation of the polynomial whose coefficients start at a(aoff).
   ! The C spells this as the POLYNOMIAL_VALUE1 macro.
   pure function poly_value1(a, aoff, order, x) result(p)
      real(dp), intent(in) :: a(0:)
      integer,  intent(in) :: aoff, order
      real(dp), intent(in) :: x
      real(dp) :: p
      integer  :: i

      p = a(aoff + order)
      do i = 1, order
         p = p * x + a(aoff + order - i)
      end do
   end function poly_value1

   ! Polish `order` roots in place, each bracketed between the previous root
   ! and the next.  Every root of a Rys polynomial is real and they interlace,
   ! so a sign change is guaranteed between consecutive brackets -- which is
   ! why failing to find one is reported rather than retried.
   function r_dnode(a, aoff, roots, order) result(err)
      real(dp), intent(in)    :: a(0:)
      integer,  intent(in)    :: aoff, order
      real(dp), intent(inout) :: roots(0:)
      integer :: err

      real(dp), parameter :: accrt = 1.0e-15_dp
      real(dp) :: x0, x1, xi, x1init, p0, p1, pi, p1init
      integer  :: m, n

      err = 0
      x1init = 0.0_dp
      p1init = a(aoff)

      do m = 0, order - 1
         x0 = x1init
         p0 = p1init
         x1init = roots(m)
         p1init = poly_value1(a, aoff, order, x1init)

         ! All coefficients zero: leave the lower-order roots alone.
         if (p1init == 0.0_dp) cycle

         if (p0 * p1init > 0.0_dp) then
            write(error_unit,'(A,I0,A,I0)') &
               "cint: root number ", m, " not found for polynomial of order ", order
            err = 1
            return
         end if

         if (x0 <= x1init) then
            x1 = x1init
            p1 = p1init
         else
            x1 = x0
            p1 = p0
            x0 = x1init
            p0 = p1init
         end if

         ! interpolate/extrapolate between [x0, x1]
         if (p1 == 0.0_dp) then
            roots(m) = x1
            cycle
         else if (p0 == 0.0_dp) then
            roots(m) = x0
            cycle
         else
            xi = x0 + (x0 - x1) / (p1 - p0) * p0
         end if

         n = 0
         do while (abs(x1 - x0) > x1 * accrt)
            n = n + 1
            if (n > 200) then
               write(error_unit,'(A)') "cint: no convergence in r_dnode"
               err = 1
               return
            end if

            pi = poly_value1(a, aoff, order, xi)
            if (pi == 0.0_dp) exit
            if (p0 * pi <= 0.0_dp) then
               x1 = xi
               p1 = pi
               xi = x0 * 0.25_dp + xi * 0.75_dp
            else
               x0 = xi
               p0 = pi
               xi = xi * 0.75_dp + x1 * 0.25_dp
            end if

            pi = poly_value1(a, aoff, order, xi)
            if (pi == 0.0_dp) exit
            if (p0 * pi <= 0.0_dp) then
               x1 = xi
               p1 = pi
            else
               x0 = xi
               p0 = pi
            end if

            xi = x0 + (x0 - x1) / (p1 - p0) * p0
         end do

         roots(m) = xi
      end do
   end function r_dnode

   ! One implicitly shifted QR sweep over the active block [n0, n1) of the
   ! Hessenberg matrix A, applied as a chain of Givens rotations.
   subroutine qr_step(a, nroots, n0, n1, shift)
      real(dp), intent(inout) :: a(0:)
      integer,  intent(in)    :: nroots, n0, n1
      real(dp), intent(in)    :: shift

      integer  :: m1, j, k, m3, j1, j2
      real(dp) :: c, s, v, x, y

      m1 = n0 + 1
      c = a(n0*nroots + n0) - shift
      s = a(m1*nroots + n0)
      v = sqrt(c*c + s*s)

      if (v == 0.0_dp) then
         v = 1.0_dp
         c = 1.0_dp
         s = 0.0_dp
      end if
      v = 1.0_dp / v
      c = c * v
      s = s * v

      do k = n0, nroots - 1
         ! Givens rotation from the left
         x = a(n0*nroots + k)
         y = a(m1*nroots + k)
         a(n0*nroots + k) = c * x + s * y
         a(m1*nroots + k) = c * y - s * x
      end do

      m3 = min(n1, n0 + 3)
      do k = 0, m3 - 1
         ! Givens rotation from the right
         x = a(k*nroots + n0)
         y = a(k*nroots + m1)
         a(k*nroots + n0) = c * x + s * y
         a(k*nroots + m1) = c * y - s * x
      end do

      do j = n0, n1 - 3
         j1 = j + 1
         j2 = j + 2
         c = a(j1*nroots + j)
         s = a(j2*nroots + j)
         v = sqrt(c*c + s*s)
         a(j1*nroots + j) = v
         a(j2*nroots + j) = 0.0_dp

         if (v == 0.0_dp) then
            v = 1.0_dp
            c = 1.0_dp
            s = 0.0_dp
         end if
         v = 1.0_dp / v
         c = c * v
         s = s * v

         do k = j1, nroots - 1
            x = a(j1*nroots + k)
            y = a(j2*nroots + k)
            a(j1*nroots + k) = c * x + s * y
            a(j2*nroots + k) = c * y - s * x
         end do

         m3 = min(n1, j + 4)
         do k = 0, m3 - 1
            x = a(k*nroots + j1)
            y = a(k*nroots + j2)
            a(k*nroots + j1) = c * x + s * y
            a(k*nroots + j2) = c * y - s * x
         end do
      end do
   end subroutine qr_step

   ! Hessenberg QR with Wilkinson shift and deflation.  Returns non-zero if a
   ! 2x2 block turns out to have complex eigenvalues, which for a Rys
   ! polynomial means the moments were already too ill-conditioned to trust.
   function hessenberg_qr(a, nroots) result(err)
      real(dp), intent(inout) :: a(0:)
      integer,  intent(in)    :: nroots
      integer :: err

      real(dp), parameter :: eps = 1.0e-15_dp
      integer,  parameter :: maxits = 30
      integer  :: n0, n1, its, k, ic, k1, m1, m2
      real(dp) :: s, a11, a22, shift, t, aa, bb

      n0 = 0
      n1 = nroots
      its = 0

      do ic = 0, nroots * maxits - 1
         k = n0
         do while (k + 1 < n1)
            s = abs(a(k*nroots + k)) + abs(a((k+1)*nroots + k + 1))
            if (abs(a((k+1)*nroots + k)) < eps * s) exit
            k = k + 1
         end do

         k1 = k + 1
         if (k1 < n1) then
            ! deflation found at position (k+1, k)
            a(k1*nroots + k) = 0.0_dp
            n0 = k1
            its = 0
            if (n0 + 1 >= n1) then
               ! block of size at most two has converged
               n0 = 0
               n1 = k1
               if (n1 < 2) then
                  err = 0
                  return
               end if
            end if
         else
            m1 = n1 - 1
            m2 = n1 - 2
            a11 = a(m1*nroots + m1)
            a22 = a(m2*nroots + m2)
            t = a11 + a22
            s = (a11 - a22)**2
            s = s + 4.0_dp * a(m1*nroots + m2) * a(m2*nroots + m1)
            if (s > 0.0_dp) then
               s = sqrt(s)
               aa = (t + s) * 0.5_dp
               bb = (t - s) * 0.5_dp
               if (abs(a11 - aa) > abs(a11 - bb)) then
                  shift = bb
               else
                  shift = aa
               end if
            else
               if (n1 == 2) then
                  write(error_unit,'(A)') "cint: hessenberg_qr failed to find real roots"
                  err = 1
                  return
               end if
               shift = t * 0.5_dp
            end if
            its = its + 1
            call qr_step(a, nroots, n0, n1, shift)
            if (its > maxits) then
               write(error_unit,'(A,I0,A)') &
                  "cint: hessenberg_qr failed to converge after ", its, " steps"
               err = 1
               return
            end if
         end if
      end do

      write(error_unit,'(A)') "cint: hessenberg_qr failed"
      err = 1
   end function hessenberg_qr

   ! cs holds the orthogonal polynomial coefficients, row `order` starting at
   ! cs(order*(nroots+1)).  The top row is turned into a companion matrix whose
   ! eigenvalues are the roots.
   function cint_polynomial_roots(roots, cs, nroots) result(err)
      real(dp), intent(inout) :: roots(0:)
      real(dp), intent(in)    :: cs(0:)
      integer,  intent(in)    :: nroots
      integer :: err

      real(dp) :: amat(0:MXRYSROOTS*MXRYSROOTS - 1)
      real(dp) :: fac, dum
      integer  :: nroots1, i, k, order

      if (nroots == 1) then
         roots(0) = -cs(2) / cs(3)
         err = 0
         return
      else if (nroots == 2) then
         dum = sqrt(cs(2*3+1)**2 - 4.0_dp*cs(2*3+0)*cs(2*3+2))
         roots(0) = (-cs(2*3+1) - dum) / cs(2*3+2) / 2.0_dp
         roots(1) = (-cs(2*3+1) + dum) / cs(2*3+2) / 2.0_dp
         err = 0
         return
      end if

      nroots1 = nroots + 1
      fac = -1.0_dp / cs(nroots*nroots1 + nroots)
      do i = 0, nroots - 1
         amat(nroots - 1 - i) = cs(nroots*nroots1 + i) * fac
      end do
      do i = nroots, nroots*nroots - 1
         amat(i) = 0.0_dp
      end do
      do i = 0, nroots - 2
         amat((i+1)*nroots + i) = 1.0_dp
      end do

      err = hessenberg_qr(amat, nroots)
      if (err == 0) then
         do i = 0, nroots - 1
            roots(nroots - 1 - i) = amat(i*nroots + i)
         end do
      else
         ! QR gave up.  Seed from the exact quadratic and polish upward, one
         ! order at a time, using the interlacing property.
         dum = sqrt(cs(2*nroots1+1) * cs(2*nroots1+1) &
                    - 4.0_dp * cs(2*nroots1+0) * cs(2*nroots1+2))
         roots(0) = 0.5_dp * (-cs(2*nroots1+1) - dum) / cs(2*nroots1+2)
         roots(1) = 0.5_dp * (-cs(2*nroots1+1) + dum) / cs(2*nroots1+2)
         do i = 2, nroots - 1
            roots(i) = 1.0_dp
         end do
         do k = 2, nroots - 1
            order = k + 1
            err = r_dnode(cs, order * nroots1, roots, order)
            if (err /= 0) exit
         end do
      end if
   end function cint_polynomial_roots

end module cint_find_roots
