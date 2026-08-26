!
! Numerical primitives for the ECP integrals.
!
! Translated from PySCF's pyscf/lib/gto/nr_ecp.c, which is where effective
! core potentials live -- libcint proper has no ECP code, so unlike every
! other module here the upstream is PySCF rather than libcint.  The licence
! is the same (Apache 2.0) and the provenance is recorded in NOTICE.
!
! Three things live here, and nothing in them depends on a basis set or a
! shell: the radial quadrature, the scaled modified spherical Bessel
! function, and the small factorial tables the rest of the ECP code indexes.
! They are separated out because they are the pieces that can be checked
! against PySCF one value at a time, before any integral is assembled.
!
! **The C's tables are not transcribed.**  nr_ecp.c carries about 4,500 lines
! of static data: a 2047-point Gauss-Chebyshev grid, and a 400x24 plus a
! 400x8x8 table of Bessel values.  Every one of them is a cache of something
! computed a few lines further down -- `ECPgauss_chebyshev` regenerates the
! grid in seventeen lines, and `ECPsph_ine` computes what the Bessel tables
! interpolate.  Porting the computation rather than the data keeps this
! module readable and costs an interpolation fast path that can be added
! later if a profile asks for it.
!
module cint_ecp_num
   use cint_const, only: dp
   implicit none
   private

   public :: ECP_LMAX, ECP_LEVEL0, ECP_LEVEL_MAX, ECP_NRS
   public :: ECP_SIM_ZERO, ECP_EXPCUTOFF, ECP_CUTOFF
   public :: RADI_POWER, SO_TYPE_OF
   public :: AS_RINV_ORIG_ATOM, AS_ECPBAS_OFFSET, AS_NECPBAS
   public :: ecp_gauss_chebyshev
   public :: ecp_sph_ine
   public :: ecp_factorial, ecp_factorial2, ecp_binom

   ! ---------------------------------------------------------------------
   ! Slots and limits, from nr_ecp.h
   ! ---------------------------------------------------------------------

   ! Extra ecpbas slots.  An ECP shell is described by the same eight-slot
   ! row a basis shell is, with two of the otherwise-unused slots carrying
   ! the r exponent and the spin-orbit flag.  ATOM_OF, ANG_OF, NPRIM_OF,
   ! NCTR_OF, PTR_EXP and PTR_COEFF keep their usual meanings, which is why
   ! ecpbas can sit inside `bas` and be reached by an offset.
   integer, parameter :: RADI_POWER = 3
   integer, parameter :: SO_TYPE_OF = 4

   !> Highest ECP angular momentum handled
   integer, parameter :: ECP_LMAX = 5

   !> env slots that say where the ECP shells are.
   !>
   !> These are PySCF's extension to the libcint env layout, not libcint's
   !> own: slots 17 to 19 are inside the reserved region below
   !> PTR_ENV_START = 20 and libcint itself never writes them.  Adopting the
   !> same numbers means an env built for PySCF works here unchanged, which
   !> is what makes a value-for-value comparison against `mole.intor` possible.
   integer, parameter :: AS_RINV_ORIG_ATOM = 17
   integer, parameter :: AS_ECPBAS_OFFSET = 18
   integer, parameter :: AS_NECPBAS = 19

   !> A coefficient below this is treated as absent
   real(dp), parameter :: ECP_SIM_ZERO = 1.0e-50_dp
   !> exp(-x) below this is zero to double precision, x = 39 giving ~1e-17
   real(dp), parameter :: ECP_EXPCUTOFF = 39.0_dp
   !> exp(x) above this overflows, x = 460 giving ~1e200
   real(dp), parameter :: ECP_CUTOFF = 460.0_dp

   !> Radial grid sizes, 2**(level+1) - 1 points.
   integer, parameter :: ECP_LEVEL0 = 5
   integer, parameter :: ECP_LEVEL_MAX = 11

   !> Points in the finest radial grid: 2047, not 2048.
   !>
   !> The C spells these two quantities the same way and they are not the
   !> same. `1 << LEVEL_MAX` is 2048 and is the size it *allocates* for the
   !> potential; the grid it integrates on is the table
   !> `rs_gauss_chebyshev2047`, which has 2047 points because the whole
   !> refinement scheme needs n+1 to be a power of two -- a strided subset of
   !> a (2^L - 1)-point rule is exactly the (2^(L-1) - 1)-point rule, which is
   !> what makes reusing the previous level's points valid and `wtscale`
   !> exact.
   !>
   !> Generating 2048 points instead gives step = 1/2049, every abscissa
   !> differs, the nesting stops being exact, and the last weight is a NaN
   !> because 1 + xi goes negative there. It converges to the same integrals
   !> anyway, which is exactly why it is worth naming rather than open-coding.
   integer, parameter :: ECP_NRS = 2**ECP_LEVEL_MAX - 1

   ! ---------------------------------------------------------------------
   ! Small tables.  Kept as data because they are exact integers up to the
   ! point where a double stops representing them exactly, and recomputing
   ! them by repeated multiplication would not reproduce the C bit for bit.
   ! ---------------------------------------------------------------------

   real(dp), parameter :: FACTORIAL(0:23) = [ &
      1.0_dp, 1.0_dp, 2.0_dp, 6.0_dp, 24.0_dp, &
      1.2e+2_dp, 7.2e+2_dp, 5.04e+3_dp, 4.032e+4_dp, 3.6288e+5_dp, &
      3.6288e+6_dp, 3.99168e+7_dp, 4.790016e+8_dp, 6.2270208e+9_dp, 8.71782912e+10_dp, &
      1.307674368e+12_dp, 2.0922789888e+13_dp, 3.55687428096e+14_dp, &
      6.402373705728e+15_dp, 1.21645100408832e+17_dp, &
      2.43290200817664e+18_dp, 5.109094217170944e+19_dp, &
      1.1240007277776077e+21_dp, 2.5852016738884978e+22_dp]

   real(dp), parameter :: FACTORIAL2(0:39) = [ &
      1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 8.0_dp, &
      15.0_dp, 48.0_dp, 105.0_dp, 384.0_dp, 945.0_dp, &
      3840.0_dp, 10395.0_dp, 46080.0_dp, 135135.0_dp, 645120.0_dp, &
      2027025.0_dp, 10321920.0_dp, 34459425.0_dp, 185794560.0_dp, 654729075.0_dp, &
      3715891200.0_dp, 13749310575.0_dp, 81749606400.0_dp, 316234143225.0_dp, &
      1961990553600.0_dp, &
      7905853580625.0_dp, 51011754393600.0_dp, 213458046676875.0_dp, &
      1428329123020800.0_dp, 6190283353629376.0_dp, &
      42849873690624000.0_dp, 1.9189878396251069e+17_dp, &
      1.371195958099968e+18_dp, 6.3326598707628524e+18_dp, &
      4.6620662575398912e+19_dp, 2.2164309547669976e+20_dp, &
      1.6783438527143608e+21_dp, 8.2007945326378929e+21_dp, &
      6.3777066403145712e+22_dp, 3.1983098677287775e+23_dp]

contains

   !> n!, and 1 outside the tabulated range the C indexes
   pure function ecp_factorial(n) result(f)
      integer, intent(in) :: n
      real(dp) :: f
      if (n < 0 .or. n > ubound(FACTORIAL, 1)) then
         f = 1.0_dp
      else
         f = FACTORIAL(n)
      end if
   end function ecp_factorial

   !> n!!, with the C's convention that a negative argument gives 1
   !>
   !> That convention is load-bearing rather than defensive: the type-2
   !> angular factors index this with (l - 1) and rely on the l = 0 case
   !> returning 1 rather than failing.
   pure function ecp_factorial2(n) result(f)
      integer, intent(in) :: n
      real(dp) :: f
      if (n < 0 .or. n > ubound(FACTORIAL2, 1)) then
         f = 1.0_dp
      else
         f = FACTORIAL2(n)
      end if
   end function ecp_factorial2

   !> Binomial coefficient n choose m
   pure function ecp_binom(n, m) result(b)
      integer, intent(in) :: n, m
      real(dp) :: b
      if (m < 0 .or. m > n) then
         b = 0.0_dp
      else
         b = ecp_factorial(n)/(ecp_factorial(m)*ecp_factorial(n - m))
      end if
   end function ecp_binom

   !> Gauss-Chebyshev radial grid on (0, 1], mapped for the ECP integrand
   !>
   !> The map is the one from Perez-Jorda, Becke and San-Fabian: a Chebyshev
   !> rule of the second kind on [-1, 1] pushed through r = 1 - log2(1 + x),
   !> which puts points where an ECP radial integrand actually has support
   !> without needing a cutoff parameter.
   !>
   !> This is the routine the C's 2047-point table was generated from; it is
   !> called instead of carrying the table.
   pure subroutine ecp_gauss_chebyshev(rs, ws, n)
      real(dp), intent(out) :: rs(0:), ws(0:)
      integer, intent(in) :: n

      real(dp), parameter :: PI = acos(-1.0_dp)
      real(dp), parameter :: LOG2E = 1.0_dp/log(2.0_dp)
      real(dp) :: step, fac, xinc, x1, x2, x3, x4, xi
      integer :: i

      step = 1.0_dp/real(n + 1, dp)
      fac = 16.0_dp*step/3.0_dp
      xinc = PI*step
      x1 = 0.0_dp
      do i = 0, n - 1
         x1 = x1 + xinc
         x2 = sin(x1)
         x3 = sin(x1*2.0_dp)
         x4 = x2*x2
         xi = real(n - i*2 - 1, dp)*step &
              + (1.0_dp/PI)*(1.0_dp + (2.0_dp/3.0_dp)*x4)*x3
         rs(i) = 1.0_dp - log(1.0_dp + xi)*LOG2E
         ws(i) = fac*x4*x4*LOG2E/(1.0_dp + xi)
      end do
   end subroutine ecp_gauss_chebyshev

   !> Scaled modified spherical Bessel function of the first kind
   !>
   !> `out(l) = i_l(z) * exp(-z)` for l = 0 .. order, which is the
   !> combination the ECP radial integrand needs and the only one that stays
   !> finite over the whole range: i_l alone overflows for large z and
   !> underflows for small.
   !>
   !> Three regimes, and the boundaries are the C's:
   !>
   !>   z < 1e-7   the leading Taylor term, z^l/(2l+1)!!, times (1 - z)
   !>   z > 16     the terminating asymptotic sum, exact for integer l
   !>   otherwise  the ascending series, summed until it stops changing
   !>
   !> The middle branch's termination test is `next == s` rather than a
   !> tolerance. That is deliberate and is kept: it stops exactly when the
   !> addition no longer moves the double, which is both the tightest
   !> available answer and reproducible against the C.
   pure subroutine ecp_sph_ine(out, order, z)
      real(dp), intent(out) :: out(0:)
      integer, intent(in) :: order
      real(dp), intent(in) :: z

      real(dp) :: z2, ti, s, next, t0
      integer :: i, k

      if (z < 1.0e-7_dp) then
         out(0) = 1.0_dp - z
         do i = 1, order
            out(i) = out(i - 1)*z/real(i*2 + 1, dp)
         end do

      else if (z > 16.0_dp) then
         z2 = -0.5_dp/z
         do i = 0, order
            ti = 0.5_dp/z
            s = ti
            do k = 1, i
               ti = ti*z2
               s = s + ti*ecp_factorial(i + k) &
                   /(ecp_factorial(k)*ecp_factorial(i - k))
            end do
            out(i) = s
         end do

      else
         z2 = 0.5_dp*z*z
         t0 = exp(-z)
         do i = 0, order
            ti = t0
            s = ti
            k = 1
            do
               ! The C is `ti *= z2 / (k * (k*2+i*2+1))`, which divides
               ! before it multiplies. Written as ti*z2/denom the rounding
               ! falls differently and the answer moves in the last few ULP,
               ! so the grouping is kept.
               ti = ti*(z2/real(k*(k*2 + i*2 + 1), dp))
               next = s + ti
               if (next == s) exit
               s = next
               k = k + 1
            end do
            t0 = t0*(z/real(i*2 + 3, dp))   ! `t0 *= z/(i*2+3)`, grouped as the C has it
            out(i) = s
         end do
      end if
   end subroutine ecp_sph_ine

end module cint_ecp_num
