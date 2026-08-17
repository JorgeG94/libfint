!
! The symmetric tridiagonal eigenproblem, translated from src/eigh.c.
!
! WHY THIS EXISTS RATHER THAN A CALL TO LAPACK
! --------------------------------------------
! libcint's eigh.c is DSTEMR from LAPACK 3.9.0, vendored and reduced to the
! tridiagonal case its Wheeler quadrature needs.  It sits behind
! `#ifdef LAPACK_FOUND`, and libcint's CMakeLists never defines that -- so the
! shipped C library always runs the vendored copy and never links LAPACK.
!
! The port originally called a real dstemr instead, which was correct but cost
! two things.  It made LAPACK the port's only external dependency, in a library
! that otherwise needs nothing but libgfortran and libm.  And it was the one
! place the port was deliberately NOT bit-identical to the C: same algorithm,
! different implementation, so the catalogue check carried a tolerance above 5
! Rys roots.
!
! Translating the C removes both.  See cint_eigh, which chooses between this
! and a real dstemr, and doc/PERFORMANCE.md for what the choice costs.
!
! The names are LAPACK's, kept deliberately: dlarrk, dlarrc, dlasq2 and the
! rest are all documented in LAPACK's own source, and renaming them to
! something more Fortran-ish would cut the reader off from that.
!
module cint_eigh_dstemr
   use cint_const, only: dp
   implicit none
   private

   public :: dstemr_diagonalize

   ! From eigh.c.  MXRYSROOTS bounds every array here; the C sizes its
   ! scratch the same way and the Wheeler caller never exceeds it.
   integer,  parameter :: MXRYSROOTS = 32
   integer,  parameter :: MAXTRY     = 6
   integer,  parameter :: MAXRQITER  = 6
   integer,  parameter :: ITMAX      = 1000
   real(dp), parameter :: FAIL_THRESHOLD = 1.0e16_dp
   real(dp), parameter :: RTOL1  = 1.5e-8_dp
   real(dp), parameter :: RTOL2  = 7.45e-11_dp
   real(dp), parameter :: MINGAP = 0.001_dp

   ! DBL_EPSILON and DBL_MIN, named as the C names them.  epsilon() and
   ! tiny() are the same values; spelling them this way keeps the
   ! correspondence with eigh.c visible line by line.
   real(dp), parameter :: DBL_EPSILON = epsilon(1.0_dp)
   real(dp), parameter :: DBL_MIN     = tiny(1.0_dp)

contains

   ! ---- the 2x2 case, dlaev2 -------------------------------------------

   integer function dlaev2(eig, vec, diag, offd) result(info)
      real(dp), intent(out) :: eig(0:), vec(0:)
      real(dp), intent(in)  :: diag(0:), offd(0:)
      real(dp) :: a, b, c, df, cs, ct, tb, sm, tn, rt, tmp
      real(dp) :: rt1, rt2, cs1, sn1
      integer  :: sgn1, sgn2

      a = diag(0); b = offd(0); c = diag(1)
      sm = a + c
      df = a - c
      tb = b + b
      rt = sqrt(tb*tb + df*df)

      if (sm > 0.0_dp) then
         rt1 = (sm + rt) * 0.5_dp
         sgn1 = 1
         rt2 = (a*c - b*b) / rt1
      else if (sm < 0.0_dp) then
         rt1 = (sm - rt) * 0.5_dp
         sgn1 = -1
         rt2 = (a*c - b*b) / rt1
      else
         rt1 = rt * 0.5_dp
         rt2 = rt * (-0.5_dp)
         sgn1 = 1
      end if

      if (df >= 0.0_dp) then
         cs = df + rt
         sgn2 = 1
      else
         cs = df - rt
         sgn2 = -1
      end if

      if (abs(cs) > abs(tb)) then
         ct = -tb / cs
         sn1 = 1.0_dp / sqrt(ct*ct + 1.0_dp)
         cs1 = ct * sn1
      else
         if (b == 0.0_dp) then
            cs1 = 1.0_dp
            sn1 = 0.0_dp
         else
            tn = -cs / tb
            cs1 = 1.0_dp / sqrt(tn*tn + 1.0_dp)
            sn1 = tn * cs1
         end if
      end if

      if (sgn1 == sgn2) then
         tmp = cs1
         cs1 = -sn1
         sn1 = tmp
      end if

      eig(0) = rt2
      eig(1) = rt1
      vec(0) = -sn1
      vec(1) = cs1
      vec(2) = cs1
      vec(3) = sn1
      info = 0
   end function dlaev2

   ! ---- bisection for one eigenvalue, dlarrk ---------------------------

   integer function dlarrk(n, iw, gl, gu, diag, e2, reltol, w, werr) result(info)
      integer,  intent(in)  :: n, iw
      real(dp), intent(in)  :: gl, gu, diag(0:), e2(0:), reltol
      real(dp), intent(out) :: w, werr
      integer  :: i, it, negcnt
      real(dp) :: mid, tmp1, left, right, tnorm

      info = 0
      if (n <= 0) return

      tnorm = max(abs(gl), abs(gu))
      info = -1
      left  = gl - tnorm * 2.0_dp * DBL_EPSILON * n
      right = gu + tnorm * 2.0_dp * DBL_EPSILON * n

      do it = 0, ITMAX - 1
         tmp1 = abs(right - left)
         if (tmp1 <= 0.0_dp .or. tmp1 < reltol * max(abs(right), abs(left))) then
            info = 0
            exit
         end if

         ! count negative pivots at the mid-point
         mid = (left + right) * 0.5_dp
         negcnt = 0
         tmp1 = diag(0) - mid
         if (tmp1 <= 0.0_dp) negcnt = negcnt + 1
         do i = 1, n - 1
            tmp1 = diag(i) - e2(i-1) / tmp1 - mid
            if (tmp1 <= 0.0_dp) negcnt = negcnt + 1
         end do
         if (negcnt >= iw) then
            right = mid
         else
            left = mid
         end if
      end do

      w    = (left + right) * 0.5_dp
      werr = abs(right - left) * 0.5_dp
   end function dlarrk

   ! ---- eigenvalue counts either side of an interval, dlarrc ----------

   subroutine dlarrc(n, vl, vu, diag, e2, lcnt, rcnt)
      integer,  intent(in)  :: n
      real(dp), intent(in)  :: vl, vu, diag(0:), e2(0:)
      integer,  intent(out) :: lcnt, rcnt
      integer  :: i
      real(dp) :: tmp, lpivot, rpivot

      lcnt = 0
      rcnt = 0
      lpivot = diag(0) - vl
      rpivot = diag(0) - vu
      if (lpivot <= 0.0_dp) lcnt = lcnt + 1
      if (rpivot <= 0.0_dp) rcnt = rcnt + 1
      do i = 0, n - 2
         tmp = e2(i)
         lpivot = diag(i+1) - vl - tmp / lpivot
         rpivot = diag(i+1) - vu - tmp / rpivot
         if (lpivot <= 0.0_dp) lcnt = lcnt + 1
         if (rpivot <= 0.0_dp) rcnt = rcnt + 1
      end do
   end subroutine dlarrc

   ! ---- the dqds shift estimate, dlasq4 --------------------------------
   !
   ! The four q/e vectors are regions of one work array and are swapped
   ! between sweeps, so they arrive as offsets into it rather than as
   ! slices.  That is not a stylistic choice: the C indexes evecp[n0-3] and
   ! friends, which for small n0 reaches BELOW the region's start into the
   ! neighbouring one.  Offsets into the whole array reproduce that exactly,
   ! where a Fortran slice would either trap or silently move the window.
   !
   ! tau is intent(inout): every early return in the C leaves *tau alone.
   subroutine dlasq4(i0, n0, n0init, w, qo, q1o, eo, e1o, dmin, dn, tau)
      integer,  intent(in)    :: i0, n0, n0init, qo, q1o, eo, e1o
      real(dp), intent(in)    :: w(0:), dmin(0:), dn(0:)
      real(dp), intent(inout) :: tau
      real(dp) :: a2, b1, b2, gap1, gap2, s
      integer  :: i

      s = 0.0_dp

      if (n0init == n0) then

         if (dmin(0) == dn(0)) then

            if (dmin(1) == dn(1)) then
               ! cases 2 and 3
               b1 = sqrt(w(qo+n0-1) * w(eo+n0-2))
               b2 = sqrt(w(qo+n0-2) * w(eo+n0-3))
               a2 = w(qo+n0-2) + w(eo+n0-2)
               gap2 = dmin(2) - a2 - dmin(2) * 0.25_dp
               if (gap2 > b2) then
                  gap1 = a2 - dn(0) - b2 / gap2 * b2
               else
                  gap1 = a2 - dn(0) - (b1 + b2)
               end if
               if (gap1 > b1) then
                  s = max(dn(0) - b1 / gap1 * b1, dmin(0) * 0.5_dp)
               else
                  s = 0.0_dp
                  if (dn(0) > b1) s = dn(0) - b1
                  if (a2 > b1 + b2) s = min(s, a2 - (b1 + b2))
                  s = max(s, dmin(0) * 0.333_dp)
               end if
            else
               ! case 4
               if (w(eo+n0-2) > w(qo+n0-2)) return
               b2 = w(eo+n0-2) / w(qo+n0-2)
               a2 = b2
               do i = n0 - 3, i0, -1
                  b1 = b2
                  if (w(eo+i) > w(qo+i)) return
                  b2 = b2 * (w(eo+i) / w(qo+i))
                  a2 = a2 + b2
                  if (0.563_dp < a2 .or. max(b1, b2) < a2 * 0.01_dp) exit
               end do
               a2 = a2 * 1.05_dp
               if (a2 < 0.563_dp) then
                  s = dn(0) * (1.0_dp - sqrt(a2)) / (a2 + 1.0_dp)
               else
                  s = dmin(0) * 0.25_dp
               end if
            end if

         else if (dmin(0) == dn(1)) then
            ! case 4
            if (w(e1o+n0-2) > w(q1o+n0-1) .or. w(eo+n0-3) > w(qo+n0-3)) return
            a2 = w(e1o+n0-2) / w(q1o+n0-1)
            b2 = w(eo+n0-3) / w(qo+n0-3)
            a2 = a2 + b2
            do i = n0 - 4, i0, -1
               if (b2 == 0.0_dp) exit
               b1 = b2
               if (w(eo+i) > w(qo+i)) return
               b2 = b2 * (w(eo+i) / w(qo+i))
               a2 = a2 + b2
               if (max(b2, b1) * 100.0_dp < a2 .or. 0.563_dp < a2) exit
            end do
            a2 = a2 * 1.05_dp
            if (a2 < 0.563_dp) then
               s = dn(1) * (1.0_dp - sqrt(a2)) / (a2 + 1.0_dp)
            else
               s = dmin(0) * 0.25_dp
            end if

         else if (dmin(0) == dn(2)) then
            ! case 5
            if (w(e1o+n0-3) > w(q1o+n0-2) .or. w(e1o+n0-2) > w(q1o+n0-1)) return
            a2 = w(e1o+n0-3) / w(q1o+n0-2) * (w(e1o+n0-2) / w(q1o+n0-1) + 1.0_dp)
            if (n0 - i0 > 3) then
               b2 = w(eo+n0-4) / w(qo+n0-4)
               a2 = a2 + b2
               do i = n0 - 5, i0, -1
                  b1 = b2
                  if (w(eo+i) > w(qo+i)) return
                  b2 = b2 * (w(eo+i) / w(qo+i))
                  a2 = a2 + b2
                  if (0.563_dp < a2 .or. max(b2, b1) < a2 * 0.01_dp) exit
               end do
               a2 = a2 * 1.05_dp
            end if
            if (a2 < 0.563_dp) then
               s = dn(2) * (1.0_dp - sqrt(a2)) / (a2 + 1.0_dp)
            else
               s = dmin(0) * 0.25_dp
            end if

         else
            ! case 6, no information
            s = dmin(0) * 0.25_dp
         end if

      else if (n0init == n0 + 1) then
         ! one eigenvalue just deflated; use dmin(1), dn(1)
         if (dmin(1) == dn(1) .and. dmin(2) == dn(2)) then
            if (w(eo+n0-2) > w(qo+n0-2)) return
            s = dmin(1) * 0.333_dp
            b1 = w(eo+n0-2) / w(qo+n0-2)
            b2 = b1
            if (b2 /= 0.0_dp) then
               do i = n0 - 3, i0, -1
                  a2 = b1
                  if (w(eo+i) > w(qo+i)) return
                  b1 = b1 * (w(eo+i) / w(qo+i))
                  b2 = b2 + b1
                  if (max(b1, a2) * 100.0_dp < b2) exit
               end do
            end if
            a2 = dmin(1) / (b2 * 1.05_dp + 1.0_dp)
            b2 = sqrt(b2 * 1.05_dp)
            gap2 = dmin(2) * 0.5_dp - a2
            if (gap2 > 0.0_dp .and. gap2 > b2 * a2) then
               s = max(s, a2 * (1.0_dp - a2 * 1.01_dp * (b2 / gap2) * b2))
            else
               s = max(s, a2 * (1.0_dp - b2 * 1.01_dp))
            end if
         else if (dmin(1) == dn(1)) then
            s = dmin(1) * 0.5_dp
         else
            s = dmin(1) * 0.25_dp
         end if

      else if (n0init == n0 + 2) then
         ! two eigenvalues deflated; use dmin(2), dn(2)
         if (dmin(2) == dn(2) .and. w(eo+n0-2) * 2.0_dp < w(qo+n0-2)) then
            if (w(eo+n0-2) > w(qo+n0-2)) return
            b1 = w(eo+n0-2) / w(qo+n0-2)
            b2 = b1
            if (b2 /= 0.0_dp) then
               do i = n0 - 2, i0 + 1, -1
                  if (w(eo+i-1) > w(qo+i-1)) return
                  b1 = b1 * (w(eo+i-1) / w(qo+i-1))
                  b2 = b2 + b1
                  if (b1 * 100.0_dp < b2) exit
               end do
            end if
            s = dmin(2) * 0.333_dp
            a2 = dmin(2) / (b2 * 1.05_dp + 1.0_dp)
            b2 = sqrt(b2 * 1.05_dp)
            gap2 = w(qo+n0-2) + w(eo+n0-3) - sqrt(w(qo+n0-3) * w(eo+n0-3)) - a2
            if (gap2 > 0.0_dp .and. gap2 > b2 * a2) then
               s = max(s, a2 * (1.0_dp - a2 * 1.01_dp * (b2 / gap2) * b2))
            else
               s = max(s, a2 * (1.0_dp - b2 * 1.01_dp))
            end if
         else
            s = dmin(2) * 0.25_dp
         end if

      else if (n0init > n0 + 2) then
         ! case 12, more than two deflated, no information
         s = 0.0_dp
      end if

      tau = s
   end subroutine dlasq4

   ! ---- one dqds sweep, dlasq5 -----------------------------------------

   subroutine dlasq5(i0, n0, w, qo, q1o, eo, e1o, tau, tol, dmin, dn)
      integer,  intent(in)    :: i0, n0, qo, q1o, eo, e1o
      real(dp), intent(inout) :: w(0:)
      real(dp), intent(in)    :: tau, tol
      real(dp), intent(out)   :: dmin(0:), dn(0:)
      real(dp) :: dg, dgmin, temp
      integer  :: j

      dg = w(qo+i0) - tau
      dgmin = dg

      do j = i0, n0 - 4
         w(q1o+j) = dg + w(eo+j)
         temp = w(qo+j+1) / w(q1o+j)
         dg = dg * temp - tau
         if (dg < tol) dg = 0.0_dp
         dgmin = min(dgmin, dg)
         w(e1o+j) = w(eo+j) * temp
      end do
      dn(2) = dg

      j = n0 - 3
      w(q1o+j) = dg + w(eo+j)
      temp = w(qo+j+1) / w(q1o+j)
      w(e1o+j) = w(eo+j) * temp
      dg = dg * temp - tau
      dn(1) = dg

      j = n0 - 2
      w(q1o+j) = dg + w(eo+j)
      temp = w(qo+j+1) / w(q1o+j)
      w(e1o+j) = w(eo+j) * temp
      dg = dg * temp - tau
      dn(0) = dg

      w(q1o+n0-1) = dg

      dmin(2) = dgmin
      dmin(1) = min(dmin(2), dn(1))
      dmin(0) = min(dmin(1), dn(0))
   end subroutine dlasq5

   ! ---- the dqds driver, dlasq2 ----------------------------------------
   !
   ! wk(wo:) is the C's `work`, four n-long regions: qvec, qvec1, evec,
   ! evec1.  The last two of each pair swap every sweep, which is why the
   ! four live as integer offsets that get exchanged rather than as arrays.
   ! qvec's FIXED offset is used as well -- the deflation writes converged
   ! eigenvalues back to qvec while reading through the swapped qvecp --
   ! so both are carried.
   integer function dlasq2(n, wk, wo, diag, offd) result(info)
      integer,  intent(in)    :: n, wo
      real(dp), intent(inout) :: wk(0:), diag(0:), offd(0:)
      integer  :: i, j, itry, iwhilb, iter
      integer  :: i0, n0, n1, n2, n0init, nbig
      integer  :: qo, q1o, eo, e1o, qop, q1op, eop, e1op, swp
      real(dp) :: emax, qmin, temp, diag_sum, tol, tol2, s, t
      real(dp) :: dmin(0:2), dn(0:2), sigma, tau

      info = 0
      dmin = 0.0_dp
      dn   = 0.0_dp

      qo  = wo
      q1o = wo + n
      eo  = wo + n*2
      e1o = wo + n*3

      j = n - 1
      do i = 0, n - 2
         temp = abs(diag(i))
         wk(q1o+i) = 0.0_dp
         wk(e1o+i) = 0.0_dp
         wk(eo+j-1) = offd(i) * offd(i) * temp
         wk(qo+j)   = temp
         j = j - 1
      end do
      wk(qo+0)    = abs(diag(n-1))
      wk(q1o+n-1) = 0.0_dp
      wk(eo+n-1)  = 0.0_dp
      wk(e1o+n-1) = 0.0_dp

      ! reverse the qd array if warranted
      if (wk(qo+0) < wk(qo+n-1) * 1.5_dp) then
         j = n - 1
         do i = 0, n/2 - 1
            temp = wk(qo+i);   wk(qo+i)   = wk(qo+j);   wk(qo+j)   = temp
            temp = wk(eo+i);   wk(eo+i)   = wk(eo+j-1); wk(eo+j-1) = temp
            j = j - 1
         end do
      end if

      ! dqd maps Z to ZZ plus Li's test
      diag_sum = wk(qo+0)
      do i = 0, n - 2
         temp = diag_sum + wk(eo+i)
         wk(q1o+i) = temp
         wk(e1o+i) = wk(qo+i+1) * (wk(eo+i) / temp)
         diag_sum  = wk(qo+i+1) * (diag_sum / temp)
      end do
      wk(q1o+n-1) = diag_sum
      diag_sum = wk(q1o+0)
      do i = 0, n - 2
         temp = diag_sum + wk(e1o+i)
         wk(qo+i) = temp
         wk(eo+i) = wk(q1o+i+1) * (wk(e1o+i) / temp)
         diag_sum = wk(q1o+i+1) * (diag_sum / temp)
      end do
      wk(qo+n-1) = diag_sum

      n0 = n
      tau = 0.0_dp
      itry = 0

      do while (n0 > 0)
         if (itry >= n) then
            info = 3
            return
         end if
         itry = itry + 1

         ! evec(n0-1) holds -sigma where the submatrix i0:n0 split off
         sigma = -wk(eo+n0-1)
         if (sigma < 0.0_dp) then
            info = 1
            return
         end if

         emax = 0.0_dp
         qmin = wk(qo+n0-1)
         i = 0
         do i = n0 - 1, 1, -1
            if (wk(eo+i-1) <= 0.0_dp) exit
            if (qmin >= emax * 4.0_dp) then
               qmin = min(qmin, wk(qo+i))
               emax = max(emax, wk(eo+i-1))
            end if
         end do
         i0 = max(i, 0)

         qop  = qo
         q1op = q1o
         eop  = eo
         e1op = e1o

         if (qmin < 0.0_dp .or. emax < 0.0_dp) then
            info = 1
            return
         end if

         dmin(0) = -max(0.0_dp, qmin - 2.0_dp * sqrt(qmin * emax))

         nbig = (n0 - i0) * 10
         iwhilb = 0
         do iwhilb = 0, nbig - 1
            n0init = n0
            tol  = DBL_EPSILON * 100.0_dp
            tol2 = tol * tol

            do while (n0 > i0)
               n1 = n0 - 1
               n2 = n0 - 2
               if (n1 == i0 .or. wk(eop+n2) < tol2 * (sigma + wk(qop+n1)) .or. &
                   wk(e1op+n2) < tol2 * wk(qop+n2)) then
                  wk(qo+n1) = wk(qop+n1) + sigma
                  n0 = n0 - 1
                  cycle
               end if
               if (n2 == i0 .or. wk(eop+n1-2) < tol2 * sigma .or. &
                   wk(e1op+n1-2) < tol2 * wk(qop+n1-2)) then
                  if (wk(qop+n1) > wk(qop+n2)) then
                     s = wk(qop+n1)
                     wk(qop+n1) = wk(qop+n2)
                     wk(qop+n2) = s
                  end if
                  t = (wk(qop+n2) - wk(qop+n1) + wk(eop+n2)) * 0.5_dp
                  if (wk(eop+n2) > wk(qop+n1) * tol2 .and. t /= 0.0_dp) then
                     s = wk(qop+n1) * (wk(eop+n2) / t)
                     s = wk(qop+n1) * (wk(eop+n2) / (t + sqrt(t * (t + s))))
                     t = wk(qop+n2) + (s + wk(eop+n2))
                     wk(qop+n1) = wk(qop+n1) * (wk(qop+n2) / t)
                     wk(qop+n2) = t
                  end if
                  wk(qo+n2) = wk(qop+n2) + sigma
                  wk(qo+n1) = wk(qop+n1) + sigma
                  n0 = n0 - 2
                  cycle
               end if
               exit
            end do
            if (n0 <= i0) exit

            if (dmin(0) <= 0.0_dp) then
               tau = -dmin(0)
            else
               call dlasq4(i0, n0, n0init, wk, qop, q1op, eop, e1op, dmin, dn, tau)
            end if

            ! call dqds until dmin > 0
            tol = DBL_EPSILON * sigma
            do iter = 0, 2
               call dlasq5(i0, n0, wk, qop, q1op, eop, e1op, tau, tol, dmin, dn)
               if (dmin(0) >= 0.0_dp) then
                  exit
               else if (dmin(1) > 0.0_dp) then
                  tau = tau + dmin(0)
               else
                  tau = tau * 0.25_dp
               end if
            end do
            if (dmin(0) < 0.0_dp) then
               tau = 0.0_dp
               call dlasq5(i0, n0, wk, qop, q1op, eop, e1op, tau, tol, dmin, dn)
            end if

            sigma = sigma + tau

            swp = qop;  qop  = q1op;  q1op = swp
            swp = eop;  eop  = e1op;  e1op = swp
         end do

         if (iwhilb >= nbig) then
            info = 2
            return
         end if
      end do
   end function dlasq2

   ! ---- the eigenvalues, via a base RRR and dqds ------------------------
   !
   ! w, werr and wgap are regions of the caller's work array, passed as
   ! offsets because the C writes wgap[-1] -- one before the region's start,
   ! which is a real element of the array before it.  wko is where this
   ! routine's own scratch begins.
   integer function compute_eigenvalues(n, diag, offd, w, wk, ero, gpo, wko) result(info)
      integer,  intent(in)    :: n, ero, gpo, wko
      real(dp), intent(inout) :: diag(0:), offd(0:), w(0:), wk(0:)
      real(dp) :: gl, gu, eabs, eold, tmp, tmp1, dmax
      real(dp) :: eps, rtol, rtl, sigma, tau, sgndef, spectral_diameter
      real(dp) :: isleft, isrght, dpivot
      integer  :: idum, ip, i, lcnt, rcnt, iinfo, e2o, wo
      logical  :: norep

      info = 0
      if (n <= 0) return

      if (n == 1) then
         w(0) = diag(0)
         wk(ero+0) = 0.0_dp
         wk(gpo+0) = 0.0_dp
         offd(0) = 0.0_dp    ! the shift for the initial RRR, zero here
         return
      end if

      eps  = 2.0_dp * DBL_EPSILON
      rtl  = 2.1e-8_dp                  ! sqrt(eps)
      rtol = 16.0_dp * DBL_EPSILON

      ! Gerschgorin interval and spectral diameter
      gl = diag(0)
      gu = diag(0)
      eold = 0.0_dp
      offd(n-1) = 0.0_dp
      do i = 0, n - 1
         wk(ero+i) = 0.0_dp
         wk(gpo+i) = 0.0_dp
         eabs = abs(offd(i))
         tmp1 = eabs + eold
         gl = min(gl, diag(i) - tmp1)
         gu = max(gu, diag(i) + tmp1)
         eold = eabs
      end do
      spectral_diameter = gu - gl
      offd(n-1) = 0.0_dp

      e2o = wko
      wo  = wko + n
      do i = 0, n - 2
         wk(e2o+i) = offd(i) * offd(i)
      end do

      ! extremal eigenvalue approximations, for the dqds case
      iinfo = dlarrk(n, 1, gl, gu, diag, wk(e2o:), rtl, tmp, tmp1)
      if (iinfo /= 0) then
         info = -1
         return
      end if
      isleft = max(gl, tmp - tmp1 - eps * 100.0_dp * abs(tmp - tmp1))
      iinfo = dlarrk(n, n, gl, gu, diag, wk(e2o:), rtl, tmp, tmp1)
      if (iinfo /= 0) then
         info = -1
         return
      end if
      isrght = min(gu, tmp + tmp1 + eps * 100.0_dp * abs(tmp + tmp1))
      spectral_diameter = isrght - isleft

      ! shift to whichever end is more populated
      call dlarrc(n, isleft + spectral_diameter * 0.25_dp, &
                     isrght - spectral_diameter * 0.25_dp, &
                  diag, wk(e2o:), lcnt, rcnt)

      if (lcnt - 1 >= n - rcnt) then
         sigma = max(isleft, gl)
         sgndef = 1.0_dp
      else
         sigma = min(isrght, gu)
         sgndef = -1.0_dp
      end if

      tau = max(spectral_diameter * n, 2.0_dp * abs(sigma)) * eps

      do idum = 0, MAXTRY - 1
         ! L D L^T of T - sigma I: D in wk(wo:), L in wk(wo+n:), reciprocal
         ! pivots in wk(wo+2n:)
         dpivot = diag(0) - sigma
         wk(wo+0) = dpivot
         dmax = abs(wk(wo+0))
         do i = 0, n - 2
            wk(wo+n*2+i) = 1.0_dp / wk(wo+i)
            tmp = offd(i) * wk(wo+n*2+i)
            wk(wo+n+i) = tmp
            dpivot = diag(i+1) - sigma - tmp * offd(i)
            wk(wo+i+1) = dpivot
            dmax = max(dmax, abs(dpivot))
         end do

         norep = dmax > spectral_diameter * 64.0_dp
         if (.not. norep) then
            ! every entry of D must share a sign
            do i = 0, n - 1
               tmp = sgndef * wk(wo+i)
               if (tmp < 0.0_dp) then
                  norep = .true.
                  exit
               end if
            end do
         end if

         if (norep) then
            if (idum == MAXTRY - 1) then
               if (sgndef == 1.0_dp) then
                  sigma = gl - spectral_diameter * 2.0_dp * eps * n
               else
                  sigma = gu + spectral_diameter * 2.0_dp * eps * n
               end if
            else
               sigma = sigma - sgndef * tau
               tau = tau * 2.0_dp
            end if
         else
            exit
         end if
      end do

      ! store the shift, then D and L
      offd(n-1) = sigma
      do ip = 0, n - 1
         diag(ip) = wk(wo+ip)
      end do
      do ip = 0, n - 2
         offd(ip) = wk(wo+n+ip)
      end do

      iinfo = dlasq2(n, wk, wo, diag, offd)
      if (iinfo /= 0) then
         info = -5
         return
      end if
      do i = 0, n - 1
         if (wk(wo+i) < 0.0_dp) then
            info = -6
            return
         end if
      end do

      if (sgndef > 0.0_dp) then
         do i = 0, n - 1
            w(i) = wk(wo+n-1-i)
         end do
      else
         do i = 0, n - 1
            w(i) = -wk(wo+i)
         end do
      end if

      do i = 0, n - 1
         wk(ero+i) = rtol * abs(w(i))
      end do

      do i = 0, n - 2
         wk(gpo+i) = max(0.0_dp, w(i+1) - wk(ero+i+1) - (w(i) + wk(ero+i)))
      end do
      ! wgap[-1] in the C: one element before the region, which is a real
      ! element of the array it sits in
      wk(gpo-1)   = max(0.0_dp, w(0) - wk(ero+0) - gl)
      wk(gpo+n-1) = max(0.0_dp, gu - sigma - (w(n-1) + wk(ero+n-1)))
   end function compute_eigenvalues

   ! ---- a shifted RRR for a cluster, dlarrf ----------------------------

   integer function dlarrf(n, diag, offd, ld, clstrt, w, wgap, werr, clgapl, &
                           sigma, dplus, lplus) result(info)
      integer,  intent(in)    :: n, clstrt
      real(dp), intent(in)    :: diag(0:), offd(0:), ld(0:), w(0:), wgap(0:), werr(0:)
      real(dp), intent(in)    :: clgapl
      real(dp), intent(out)   :: sigma, dplus(0:), lplus(0:)
      integer  :: i, ktry
      real(dp) :: s, tmp, max1, growthbound, lsigma

      info = 0
      ! a small fudge, to be sure of shifting to the outside
      lsigma = w(clstrt) - werr(clstrt)
      lsigma = lsigma - abs(lsigma) * 4.0_dp * DBL_EPSILON
      growthbound = diag(0) * 8.0_dp
      max1 = 0.0_dp

      do ktry = 0, 1
         s = -lsigma
         dplus(0) = diag(0) + s
         max1 = abs(dplus(0))
         do i = 0, n - 2
            tmp = ld(i) / dplus(i)
            lplus(i) = tmp
            s = s * tmp * offd(i) - lsigma
            dplus(i+1) = diag(i+1) + s
            max1 = max(max1, abs(dplus(i+1)))
         end do
         sigma = lsigma
         if (max1 <= growthbound) return
         ! failed the RRR test; back off to the outside
         lsigma = lsigma - min(clgapl * 0.25_dp, wgap(clstrt))
      end do

      if (max1 > FAIL_THRESHOLD) info = 1
   end function dlarrf

   ! ---- negative-pivot count at a shift, dlaneg -----------------------

   integer function dlaneg(n, diag, lld, sigma, twist_index) result(negcnt)
      integer,  intent(in) :: n, twist_index
      real(dp), intent(in) :: diag(0:), lld(0:), sigma
      integer  :: j
      real(dp) :: p, t, dplus, dminus

      negcnt = 0
      ! upper part: L D L^T - sigma I = L+ D+ L+^T
      t = -sigma
      do j = 0, twist_index - 2
         dplus = diag(j) + t
         if (dplus < 0.0_dp) negcnt = negcnt + 1
         t = t / dplus * lld(j) - sigma
      end do

      ! lower part: L D L^T - sigma I = U- D- U-^T
      p = diag(n-1) - sigma
      do j = n - 2, twist_index - 1, -1
         dminus = lld(j) + p
         if (dminus < 0.0_dp) negcnt = negcnt + 1
         p = p / dminus * diag(j) - sigma
      end do

      ! the twist index; T was shifted by sigma initially
      if (t + sigma + p < 0.0_dp) negcnt = negcnt + 1
   end function dlaneg

   ! ---- bisection refinement, dlarrb -----------------------------------

   integer function dlarrb(n, diag, lld, ifirst, ilast, rtol1_, rtol2_, &
                           w, wgap, werr, twist_index) result(info)
      integer,  intent(in)    :: n, ifirst, ilast, twist_index
      real(dp), intent(in)    :: diag(0:), lld(0:), rtol1_, rtol2_
      real(dp), intent(inout) :: w(0:), wgap(0:), werr(0:)
      integer  :: i, iter, negcnt
      real(dp) :: mid, back, left, right, width, cvrgd

      info = 0
      mid = 0.0_dp
      width = 0.0_dp
      cvrgd = max(rtol1_ * wgap(ifirst), &
                  rtol2_ * max(abs(w(ifirst)), abs(w(ilast-1))))

      do i = ifirst, ilast - 1
         if (werr(i) < cvrgd) cycle

         ! bracket the eigenvalue: negcount from L+D+L+^T = L D L^T - left
         left = w(i)
         back = werr(i)
         negcnt = ilast
         do while (negcnt > i)
            left = left - back
            back = back * 2.0_dp
            negcnt = dlaneg(n, diag, lld, left, twist_index)
         end do

         right = w(i)
         back = werr(i)
         negcnt = ifirst
         do while (negcnt <= i)
            right = right + back
            back = back * 2.0_dp
            negcnt = dlaneg(n, diag, lld, right, twist_index)
         end do

         do iter = 0, ITMAX - 1
            mid = (left + right) * 0.5_dp
            width = right - mid
            if (width < cvrgd) exit
            negcnt = dlaneg(n, diag, lld, mid, twist_index)
            if (negcnt <= i) then
               left = mid
            else
               right = mid
            end if
         end do
         w(i) = mid
         werr(i) = width
      end do

      do i = ifirst, ilast - 2
         wgap(i) = max(0.0_dp, w(i+1) - werr(i+1) - w(i) - werr(i))
      end do
   end function dlarrb

   ! ---- one eigenvector by the twisted factorisation, dlar1v ----------

   subroutine dlar1v(n, lambda, diag, offd, ld, lld, gaptol, vec, negcnt, &
                     twist_index, resid, rqcorr, wk, wo)
      integer,  intent(in)    :: n, wo
      real(dp), intent(in)    :: lambda, diag(0:), offd(0:), ld(0:), lld(0:), gaptol
      real(dp), intent(inout) :: vec(0:), wk(0:)
      integer,  intent(out)   :: negcnt
      integer,  intent(inout) :: twist_index
      real(dp), intent(out)   :: resid, rqcorr
      integer  :: i, r1, r2, neg1, neg2, lpo, umo, wpo
      real(dp) :: s, tmp, nrminv, mingma, dplus, dminus, ztz

      lpo = wo            ! lplus
      umo = wo + n        ! uminus
      wpo = wo + n*2      ! work_p

      if (twist_index == -1) then
         r1 = 0
         r2 = n
         twist_index = 0
      else
         r1 = twist_index
         r2 = twist_index + 1
      end if

      ! progressive transform, differential form, down to r1
      neg2 = 0
      s = diag(n-1) - lambda
      wk(wpo+n-1) = s
      do i = n - 2, r1, -1
         dminus = lld(i) + s
         if (dminus < 0.0_dp) neg2 = neg2 + 1
         tmp = diag(i) / dminus
         wk(umo+i) = offd(i) * tmp
         s = s * tmp - lambda
         wk(wpo+i) = s
      end do

      ! stationary transform, differential form, up to r2
      neg1 = 0
      s = -lambda
      do i = 0, r1 - 1
         dplus = diag(i) + s
         if (dplus < 0.0_dp) neg1 = neg1 + 1
         wk(lpo+i) = ld(i) / dplus
         s = s * wk(lpo+i) * offd(i) - lambda
      end do
      mingma = s + lambda + wk(wpo+r1)
      if (mingma < 0.0_dp) neg1 = neg1 + 1

      negcnt = neg1 + neg2

      ! largest diagonal element of the inverse, between r1 and r2
      do i = r1, r2 - 2
         dplus = diag(i) + s
         wk(lpo+i) = ld(i) / dplus
         tmp = s * wk(lpo+i) * offd(i)
         s = tmp - lambda
         tmp = tmp + wk(wpo+i+1)
         if (abs(tmp) <= abs(mingma)) then
            mingma = tmp
            twist_index = i + 1
         end if
      end do

      ! the FP vector: solve N^T v = e_r
      vec(twist_index) = 1.0_dp
      ztz = 1.0_dp
      do i = twist_index - 1, 0, -1
         tmp = -(wk(lpo+i) * vec(i+1))
         ztz = ztz + tmp * tmp
         vec(i) = tmp
      end do
      do i = twist_index, n - 2
         tmp = -(wk(umo+i) * vec(i))
         ztz = ztz + tmp * tmp
         vec(i+1) = tmp
      end do

      tmp = 1.0_dp / ztz
      nrminv = sqrt(tmp)
      do i = 0, n - 1
         vec(i) = vec(i) * nrminv
      end do
      resid  = abs(mingma) * nrminv
      rqcorr = mingma * tmp
   end subroutine dlar1v

   ! ---- the eigenvectors, by MRRR --------------------------------------
   !
   ! wgap is reached through an offset because the C reads wgap[-1] when a
   ! cluster starts at index 0 -- see the note on compute_eigenvalues.
   ! vec doubles as scratch: a cluster's RRR is parked in the rows of vec
   ! belonging to its leftmost eigenvalue until that cluster is resolved.
   integer function compute_eigenvectors(n, diag, offd, w, wk, ero, gpo, wko, &
                                         vec, iwk) result(info)
      integer,  intent(in)    :: n, ero, gpo, wko
      real(dp), intent(inout) :: diag(0:), offd(0:), w(0:), wk(0:), vec(0:)
      integer,  intent(inout) :: iwk(0:)
      integer  :: i, j, k, icls, iter, idone, ndepth
      integer  :: ncluster, ncluster1, negcnt
      integer  :: oldfst, oldlst, newfst, newlst, iinfo
      integer  :: bwo, ldo, lldo, wrko, tio, oclo, nclo, swp
      logical  :: needbs
      real(dp) :: fudge, eps, rqtol, tol, tmp
      real(dp) :: left, right, gap, bstw, savgap, gaptol
      real(dp) :: sigma, tau, resid, lambda, bstres
      real(dp) :: rqcorr, resid_tol, rqcorr_tol

      info = 0
      if (n <= 0) return

      ! the first n entries of the scratch hold the shifted eigenvalues
      bwo  = wko
      ldo  = wko + n
      lldo = wko + n*2
      wrko = wko + n*3
      do i = 0, n*6 - 1
         wk(wko+i) = 0.0_dp
      end do
      tio  = 0            ! twist indices
      oclo = n            ! old cluster ranges
      nclo = n*3          ! new cluster ranges
      do i = 0, n - 1
         iwk(tio+i) = 0
      end do

      eps   = DBL_EPSILON
      rqtol = DBL_EPSILON * 2.0_dp
      tol   = DBL_EPSILON * 8.0_dp

      sigma = offd(n-1)
      if (n == 1) then
         vec(0) = 1.0_dp
         w(0) = w(0) + sigma
         return
      end if

      do i = 0, n - 1
         wk(bwo+i) = w(i)
         w(i) = w(i) + sigma
      end do

      ncluster = 1
      iwk(oclo+0) = 0
      iwk(oclo+1) = n
      idone = 0

      do ndepth = 0, n - 1
         if (idone == n) exit

         ncluster1 = ncluster
         ncluster = 0
         do icls = 0, ncluster1 - 1
            oldfst = iwk(oclo + icls*2)
            oldlst = iwk(oclo + icls*2 + 1)
            if (ndepth > 0) then
               ! the RRR computed at the previous level, parked in vec
               do i = 0, n - 1
                  diag(i) = vec(oldfst*n + i)
               end do
               do i = 0, n - 2
                  offd(i) = vec((oldfst+1)*n + i)
               end do
               sigma = vec((oldfst+2)*n - 1)
            end if

            ! DL and DLL of the current RRR
            do j = 0, n - 2
               tmp = diag(j) * offd(j)
               wk(ldo+j)  = tmp
               wk(lldo+j) = tmp * offd(j)
            end do

            newfst = oldfst
            do newlst = oldfst + 1, oldlst
               if (newlst < oldlst) then
                  if (wk(gpo+newlst-1) < MINGAP * abs(wk(bwo+newlst-1))) cycle
               end if

               if (newlst - newfst > 1) then
                  ! a cluster: shift as close as possible to get large
                  ! relative gaps, then hand the child down a level
                  iinfo = dlarrb(n, diag, wk(lldo:), newfst, newfst+1, rqtol, rqtol, &
                                 wk(bwo:), wk(gpo:), wk(ero:), n)
                  iinfo = dlarrb(n, diag, wk(lldo:), newlst-1, newlst, rqtol, rqtol, &
                                 wk(bwo:), wk(gpo:), wk(ero:), n)

                  iinfo = dlarrf(n, diag, offd, wk(ldo:), newfst, wk(bwo:), &
                                 wk(gpo:), wk(ero:), wk(gpo+newfst-1), &
                                 tau, vec(newfst*n:), vec((newfst+1)*n:))
                  if (iinfo /= 0) then
                     info = -2
                     return
                  end if
                  vec((newfst+2)*n - 1) = sigma + tau

                  do k = newfst, newlst - 1
                     fudge = eps * 3.0_dp * abs(wk(bwo+k))
                     wk(bwo+k) = wk(bwo+k) - tau
                     fudge = fudge + eps * 4.0_dp * abs(wk(bwo+k))
                     wk(ero+k) = wk(ero+k) + fudge
                  end do
                  iwk(nclo + ncluster*2)     = newfst
                  iwk(nclo + ncluster*2 + 1) = newlst
                  ncluster = ncluster + 1

               else
                  ! a singleton: Rayleigh quotient iteration, falling back
                  ! to bisection when the correction points the wrong way
                  lambda = wk(bwo+newfst)
                  left  = lambda - wk(ero+newfst)
                  right = lambda + wk(ero+newfst)
                  gap = wk(gpo+newfst)
                  if (newfst == 0 .or. newfst+1 == n) then
                     gaptol = 0.0_dp
                  else
                     gaptol = gap * eps
                  end if
                  savgap = gap
                  resid_tol  = tol * gap
                  rqcorr_tol = rqtol * abs(lambda)

                  needbs = .false.
                  bstres = 1.0e307_dp
                  bstw = 0.0_dp
                  iwk(tio+newfst) = -1
                  do iter = 0, MAXRQITER - 1
                     call dlar1v(n, lambda, diag, offd, wk(ldo:), wk(lldo:), gaptol, &
                                 vec(newfst*n:), negcnt, iwk(tio+newfst), &
                                 resid, rqcorr, wk, wrko)

                     if (resid < bstres) then
                        bstres = resid
                        bstw = lambda
                     end if

                     if (resid < resid_tol .or. abs(rqcorr) < rqcorr_tol) exit

                     ! only use the correction when it improves the iterate
                     if (lambda + rqcorr > right .or. lambda + rqcorr < left) then
                        needbs = .true.
                        exit
                     end if

                     ! and only when it does not walk towards a neighbour
                     if (newfst < negcnt) then
                        if (rqcorr > 0.0_dp) then
                           needbs = .true.
                           exit
                        end if
                        right = lambda
                     else
                        if (rqcorr < 0.0_dp) then
                           needbs = .true.
                           exit
                        end if
                        left = lambda
                     end if

                     wk(bwo+newfst) = (right + left) * 0.5_dp
                     lambda = lambda + rqcorr
                     wk(ero+newfst) = (right - left) * 0.5_dp

                     if (right - left < rqcorr_tol) then
                        if (bstres < resid) then
                           lambda = bstw
                           call dlar1v(n, lambda, diag, offd, wk(ldo:), wk(lldo:), &
                                       gaptol, vec(newfst*n:), negcnt, &
                                       iwk(tio+newfst), resid, rqcorr, wk, wrko)
                        end if
                        exit
                     end if
                  end do

                  if (needbs) then
                     iinfo = dlarrb(n, diag, wk(lldo:), newfst, newfst+1, 0.0_dp, &
                                    eps * 2.0_dp, wk(bwo:), wk(gpo:), wk(ero:), &
                                    iwk(tio+newfst) + 1)
                     lambda = wk(bwo+newfst)
                     ! reset the twist index: an inaccurate lambda has to be
                     ! made to recompute the true mingma
                     iwk(tio+newfst) = -1
                     call dlar1v(n, lambda, diag, offd, wk(ldo:), wk(lldo:), gaptol, &
                                 vec(newfst*n:), negcnt, iwk(tio+newfst), &
                                 resid, rqcorr, wk, wrko)
                  end if

                  w(newfst) = lambda + sigma

                  ! gaps may only grow: shrinking reflects cancellation
                  ! rather than the theory
                  if (newfst > 0) then
                     wk(gpo+newfst-1) = max(wk(gpo+newfst-1), &
                        w(newfst) - wk(ero+newfst) - w(newfst-1) - wk(ero+newfst-1))
                  end if
                  if (newfst < n-1) then
                     wk(gpo+newfst) = max(savgap, &
                        w(newfst+1) - wk(ero+newfst+1) - w(newfst) - wk(ero+newfst))
                  end if

                  idone = idone + 1
               end if

               newfst = newlst
            end do
         end do

         swp = oclo; oclo = nclo; nclo = swp
      end do

      if (idone < n) info = -2
   end function compute_eigenvectors

   ! ---- the entry point, _CINTdiagonalize ------------------------------
   !
   ! The scratch layout is the C's: one real work array carved into the
   ! eigenvalues, their errors, the gaps and the routines' own space, and
   ! one integer array for the twist indices and cluster ranges.  The odd
   ! +1 on two of the offsets is what gives wgap[-1] somewhere real to
   ! live.
   integer function dstemr_diagonalize(n, diag, offd, eig, vec) result(info)
      integer,  intent(in)    :: n
      real(dp), intent(inout) :: diag(0:), offd(0:)
      real(dp), intent(inout) :: eig(0:), vec(0:)
      real(dp) :: wk(0:MXRYSROOTS*9)
      integer  :: iwk(0:MXRYSROOTS*5-1)
      integer  :: ero, gpo, wko

      info = 0
      if (n == 0) return
      if (n == 1) then
         eig(0) = diag(0)
         vec(0) = 1.0_dp
         return
      end if
      if (n == 2) then
         info = dlaev2(eig, vec, diag, offd)
         return
      end if

      ero = n            ! werr
      gpo = n*2 + 1      ! wgap, offset so that wgap(-1) is in range
      wko = n*3 + 1      ! the routines' own scratch

      info = compute_eigenvalues(n, diag, offd, eig, wk, ero, gpo, wko)
      if (info == 0) then
         info = compute_eigenvectors(n, diag, offd, eig, wk, ero, gpo, wko, vec, iwk)
      end if
   end function dstemr_diagonalize

end module cint_eigh_dstemr
