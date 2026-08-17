!
! Primitive-pair screening.
!
! Ported from the parts of src/optimizer.c that the drivers need directly.
! The rest of that file -- the CINTOpt cache with its ragged pointer arrays --
! belongs to D8; these three routines are the ones CINT1e_loop calls whether or
! not an optimizer was supplied.
!
! Note `approx_log` in the C is `#define approx_log log`; the commented-out
! fast approximation above it is not used.
!
module cint_screen
   use cint_const, only: dp
   use cint_envs,  only: pair_data, cint_opt_t, OPT_NOVALUE, LMAX1, &
                         MAX_PGTO_FOR_PAIRDATA
   implicit none
   private

   public :: pair_data, cint_set_pairdata
   public :: cint_log_max_pgto_coeff, cint_non0coeff_byshell
   public :: cint_opt_t, OPT_NOVALUE, LMAX1, MAX_PGTO_FOR_PAIRDATA

contains

   ! Largest contraction coefficient on each primitive, logged.  Used to bound
   ! a pair's contribution before computing it.
   pure subroutine cint_log_max_pgto_coeff(log_maxc, coeff, coff, nprim, nctr)
      real(dp), intent(out) :: log_maxc(0:)
      real(dp), intent(in)  :: coeff(0:)
      integer,  intent(in)  :: coff, nprim, nctr
      integer  :: i, ip
      real(dp) :: maxc
      do ip = 0, nprim - 1
         maxc = 0.0_dp
         do i = 0, nctr - 1
            maxc = max(maxc, abs(coeff(coff + i*nprim + ip)))
         end do
         log_maxc(ip) = log(maxc)
      end do
   end subroutine cint_log_max_pgto_coeff

   ! For each primitive, the contraction indices with a non-zero coefficient
   ! first, then the zero ones.  prim_to_ctr_1 walks only the first non0ctr of
   ! them; prim_to_ctr_0 needs the zeros present too, which is why they are
   ! appended rather than dropped.
   pure subroutine cint_non0coeff_byshell(sortedidx, non0ctr, ci, coff, iprim, ictr)
      integer,  intent(out) :: sortedidx(0:), non0ctr(0:)
      real(dp), intent(in)  :: ci(0:)
      integer,  intent(in)  :: coff, iprim, ictr
      integer :: ip, j, k, base

      ! TWO PASSES, NO SCRATCH.  The C collects the zero-coefficient indices
      ! in a C99 VLA -- `FINT zeroidx[ictr]` -- which costs it a stack
      ! adjustment.  The Fortran equivalent is an automatic array with a
      ! runtime bound, and gfortran puts those on the heap: this routine was
      ! calling malloc and free 414,720 times in one benchmark, four per
      ! shell quartet, for about 75 million instructions.
      !
      ! Walking the coefficients twice costs a second comparison per
      ! contraction and no allocation at all.  The order is unchanged: the
      ! second pass appends the zeros in the order it meets them, which is
      ! the order the scratch array held them.
      base = 0
      do ip = 0, iprim - 1
         k = 0
         do j = 0, ictr - 1
            if (ci(coff + iprim*j + ip) /= 0.0_dp) then
               sortedidx(base + k) = j
               k = k + 1
            end if
         end do
         non0ctr(ip) = k
         do j = 0, ictr - 1
            if (ci(coff + iprim*j + ip) == 0.0_dp) then
               sortedidx(base + k) = j
               k = k + 1
            end if
         end do
         base = base + ictr
      end do
   end subroutine cint_non0coeff_byshell

   ! Centre and exponential prefactor of every primitive pair, plus a bound
   ! (cceij) on its size.  Returns .true. when every pair is negligible, which
   ! lets the caller skip the shell pair entirely.
   !
   ! The bound comes from the overlap estimate in the C's comment:
   !   (aj*d/sqrt(aij)+1)^li (ai*d/sqrt(aij)+1)^lj pi^1.5/aij^((li+lj+3)/2)
   !     * exp(-ai*aj/aij*rr_ij)
   ! loosened to something cheap to evaluate.
   function cint_set_pairdata(pd, ai, aioff, aj, ajoff, ri, rj, &
                              log_maxci, log_maxcj, li_ceil, lj_ceil, &
                              iprim, jprim, rr_ij, expcutoff, omega) result(empty)
      type(pair_data), intent(inout) :: pd(0:)
      real(dp), intent(in) :: ai(0:), aj(0:), ri(0:), rj(0:)
      real(dp), intent(in) :: log_maxci(0:), log_maxcj(0:)
      integer,  intent(in) :: aioff, ajoff, li_ceil, lj_ceil, iprim, jprim
      real(dp), intent(in) :: rr_ij, expcutoff, omega
      logical :: empty

      integer  :: ip, jp, n, lij
      real(dp) :: aij, eij, cceij, wj, log_rr_ij, dist_ij, omega2, theta
      real(dp), parameter :: R_GUESS = 8.0_dp

      aij = ai(aioff + iprim - 1) + aj(ajoff + jprim - 1)
      log_rr_ij = 1.7_dp - 1.5_dp * log(aij)
      lij = li_ceil + lj_ceil
      if (lij > 0) then
         dist_ij = sqrt(rr_ij)
         if (omega < 0.0_dp) then
            ! attenuated kernel: the effective range is shorter, so the bound
            ! may be tightened by the screening length
            omega2 = omega * omega
            theta = omega2 / (omega2 + aij)
            log_rr_ij = log_rr_ij + lij * log(dist_ij + theta*R_GUESS + 1.0_dp)
         else
            log_rr_ij = log_rr_ij + lij * log(dist_ij + 1.0_dp)
         end if
      end if

      empty = .true.
      n = 0
      do jp = 0, jprim - 1
         do ip = 0, iprim - 1
            aij = 1.0_dp / (ai(aioff+ip) + aj(ajoff+jp))
            eij = rr_ij * ai(aioff+ip) * aj(ajoff+jp) * aij
            cceij = eij - log_rr_ij - log_maxci(ip) - log_maxcj(jp)
            pd(n)%cceij = cceij
            if (cceij < expcutoff) then
               empty = .false.
               wj = aj(ajoff+jp) * aij
               pd(n)%rij(0) = ri(0) + wj * (rj(0) - ri(0))
               pd(n)%rij(1) = ri(1) + wj * (rj(1) - ri(1))
               pd(n)%rij(2) = ri(2) + wj * (rj(2) - ri(2))
               pd(n)%eij = exp(-eij)
            else
               ! far enough away to be dropped; the coordinates are poisoned
               ! rather than left stale so a bug reading them is obvious
               pd(n)%rij = 1.0e18_dp
               pd(n)%eij = 0.0_dp
            end if
            n = n + 1
         end do
      end do
   end function cint_set_pairdata

end module cint_screen
