! -*- fortran -*-
!
! Schmidt orthogonalisation of the Rys moments, and the RDK root/weight solve
! built on it.  Included once per precision by cint_schmidt.F90, into a module
! that has bound `rk` and pulled in the matching Boys function.
!
! The C writes this three times (R_dsmit / R_lsmit / R_qsmit, each with its own
! copy of the driver); two instantiations cover it here, per
! PORT_TO_FORTRAN.md 3.1.
!
! `cs` is conceptually (n x n), flat and 0-based with the C's column-major
! addressing cs(i + k*n) preserved exactly.  Note that this is the one place
! the C indexes a table column-first; keeping it verbatim is safer than
! tidying it.

   ! Gram-Schmidt on the moment sequence.  Returns 0 on success, or the order
   ! at which the norm went non-positive -- which means the moments have lost
   ! enough precision that the higher Rys polynomials are meaningless, and the
   ! caller treats their contribution as zero.
   function r_smit(cs, fmt_ints, n) result(err)
      RKTYPE, intent(inout) :: cs(0:)
      RKTYPE, intent(in)    :: fmt_ints(0:)
      integer,  intent(in)    :: n
      integer :: err

      integer  :: i, j, k
      RKTYPE :: fac, dot, tmp
      RKTYPE :: v(0:MXRYSROOTS-1)

      err = 0
      fac = -fmt_ints(1) / fmt_ints(0)
      tmp = fmt_ints(2) + fac * fmt_ints(1)
      if (tmp <= 0.0_rk) then
         call set_zero(cs, n, 1)
         err = 1
         return
      end if
      tmp = 1.0_rk / sqrt(tmp)
      cs(0 + 0*n) = 1.0_rk / sqrt(fmt_ints(0))
      cs(0 + 1*n) = fac * tmp
      cs(1 + 1*n) = tmp

      do j = 2, n - 1
         do k = 0, j - 1
            v(k) = 0.0_rk
         end do
         fac = fmt_ints(j + j)
         do k = 0, j - 1
            dot = 0.0_rk
            do i = 0, k
               dot = dot + cs(i + k*n) * fmt_ints(i + j)
            end do
            do i = 0, k
               v(i) = v(i) - dot * cs(i + k*n)
            end do
            fac = fac - dot * dot
         end do

         if (fac <= 0.0_rk) then
            call set_zero(cs, n, j)
            if (fac == 0.0_rk) then
               err = 0
            else
               err = j
            end if
            return
         end if
         fac = 1.0_rk / sqrt(fac)
         cs(j + j*n) = fac
         do k = 0, j - 1
            cs(k + j*n) = fac * v(k)
         end do
      end do
   end function r_smit

   pure subroutine set_zero(a, n, start)
      RKTYPE, intent(inout) :: a(0:)
      integer,  intent(in)    :: n, start
      integer :: i, k
      do k = start, n - 1
         do i = 0, n - 1
            a(i + k*n) = 0.0_rk
         end do
      end do
   end subroutine set_zero

   ! Roots and weights by the RDK scheme: orthogonalise the moments, take the
   ! roots of the top polynomial, then weight each by the reciprocal of the
   ! summed squared polynomial values there.
   !
   ! The orthogonalisation runs at `rk`, but the root finding always drops to
   ! double, exactly as the C does -- _CINT_polynomial_roots has no extended
   ! variant.
   function rys_schmidt(nroots, x, lower, roots, weights) result(err)
      integer,  intent(in)  :: nroots
      real(dp), intent(in)  :: x, lower
      real(dp), intent(out) :: roots(0:), weights(0:)
      integer :: err

      RKTYPE :: fmt_ints(0:2*MXRYSROOTS-1)
      RKTYPE :: qcs(0:MXRYSROOTS*MXRYSROOTS-1)
      real(dp) :: rt(0:MXRYSROOTS-1)
      real(dp) :: cs(0:MXRYSROOTS*MXRYSROOTS-1)
      real(dp) :: root, poly, dum, dum0
      integer  :: nroots1, i, j, k

      err = 0
      nroots1 = nroots + 1

      if (lower == 0.0_dp) then
         call gamma_inc_like(fmt_ints, TO_RK(x), nroots*2)
      else
         call fmt_erfc_like(fmt_ints, TO_RK(x), TO_RK(lower), nroots*2)
      end if

      ! Guard against the moments underflowing to nothing; without it the
      ! normalisation below divides by zero.
      if (fmt_ints(0) == 0.0_rk) then
         do k = 0, nroots - 1
            roots(k) = 0.0_dp
            weights(k) = 0.0_dp
         end do
         return
      end if

      cs = 0.0_dp

      ! And again here.  The double copy returns f1/(f0-f1) directly; the
      ! extended copies form f1/f0 and let the shared tail apply
      ! root/(1-root).  Algebraically the same, numerically not: the two
      ! round differently, which shows up at the 1e-15 level.
      if (nroots == 1 .and. one_root_early_return) then
         roots(0) = TO_DP(fmt_ints(1) / (fmt_ints(0) - fmt_ints(1)))
         weights(0) = TO_DP(fmt_ints(0))
         return
      end if

      if (nroots == 1) then
         rt(0) = TO_DP(fmt_ints(1) / fmt_ints(0))
      else
         qcs = 0.0_rk
         err = r_smit(qcs, fmt_ints, nroots1)
         if (err /= 0) then
            ! The C's copies disagree here: the double one normalises any
            ! r_smit failure to 1, the extended ones pass back the order at
            ! which orthogonality was lost.  Callers only test for non-zero,
            ! so it makes no difference to behaviour -- but reproducing each
            ! copy's own convention is what lets the double path be compared
            ! against the C bit for bit.
            if (smit_err_as_one) err = 1
            return
         end if
         do k = 1, nroots1 - 1
            do i = 0, k
               cs(k*nroots1 + i) = TO_DP(qcs(k*nroots1 + i))
            end do
         end do
         err = cint_polynomial_roots(rt, cs, nroots)
         if (err /= 0) return
      end if

      dum0 = TO_DP(1.0_rk / fmt_ints(0))
      do k = 0, nroots - 1
         root = rt(k)
         ! A root pinned at 1 marks a polynomial that r_smit gave up on; its
         ! contribution is taken as negligible rather than propagated as an
         ! infinity through root/(1-root).
         if (root == 1.0_dp) then
            roots(k) = 0.0_dp
            weights(k) = 0.0_dp
            cycle
         end if

         dum = dum0
         do j = 1, nroots - 1
            poly = poly_value1(cs, j*nroots1, j, root)
            dum = dum + poly * poly
         end do
         roots(k) = root / (1.0_dp - root)
         weights(k) = 1.0_dp / dum
      end do
   end function rys_schmidt
