!
! Double-double arithmetic: ~106 bits of significand out of two real64s.
!
! Why this exists.  The extended ladder in cint_const uses real128, and
! nvfortran and LLVM Flang do not have it -- `iso_fortran_env::real128` is -1
! there and `real(real128)` is a compile error, so libfint cannot be built with
! either.  The quad precision itself is not the problem: nvc's __float128 is a
! true binary128 with a 113-bit significand, measured, identical to gcc's.  The
! gap is the Fortran front end alone.
!
! A C shim cannot fix that, because what is missing is a *kind parameter*: every
! declaration, literal and intrinsic would have to become a call.  Two real64s
! carrying a value and its rounding error do fix it, in plain Fortran that every
! compiler has.
!
! **The algorithms are exact, and that is the point.**  `two_sum` and
! `two_prod` return the rounded result *and* the error it discarded, with no
! approximation: the error of a floating-point sum is itself representable.
! Dekker splitting is used for the product rather than an FMA, because Fortran
! has no portable FMA before F2018's IEEE_FMA and nvfortran's support for it is
! not something to rely on here.
!
! **The flag contract, measured rather than assumed.**  Folding the raw bits of
! 2000 dd results over sqrt/exp/erf/div/mul:
!
!   gfortran -O2                     hi 0854955B40D86FC2   lo 85F4EDB9AB7A5253
!   nvfortran -O2 -Mnofma            hi 0854955B40D86FC2   lo 85F4EDB9AB7A5253
!   nvfortran -O2   (FMA is default) hi 0854955B40D86FC2   lo B7E34D1CF97CFA3B
!   nvfortran -fast                  hi 0854955B40D86FC2   lo B7E34D1CF97CFA3B
!   gfortran -O2 -ffast-math         hi 1B7A0DE5351EEF18   lo 0000000000000000
!
! Two different failures, and they need different responses.
!
! **FMA contraction costs reproducibility, not accuracy.**  `hi` is identical
! everywhere; only the error term changes, and it stays a valid error term --
! accuracy over [-35,35] is 3.2e-30 contracted against 3.3e-30 not.  So an FMA
! build is *correct* and merely disagrees bit-for-bit with every other build.
! nvfortran contracts by default at -O2, so -Mnofma is required to match.  On
! ifx the spelling is -no-fma; -fno-fma is accepted and silently ignored.
!
! **-ffast-math destroys it outright.**  Every `lo` folds to exactly zero: the
! reassociation of `(a + b) - a` that the error terms are built from is gone,
! and what is left is double precision wearing a double-double type.  That is
! the one that must never be allowed, and dd_check catches it.
!
! Both findings are from this repository's own compilers, not from lore; the
! landmine catalogue in bitwise_adventures documents the same FMA behaviour
! from an independent direction.
!
module cint_dd
   use iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter :: dp = real64

   type, public :: dd
      real(dp) :: hi = 0.0_dp
      real(dp) :: lo = 0.0_dp
   end type dd

   public :: dd_from, dd_to_dp, two_sum, two_prod, dd_add, dd_mul
   public :: dd_sub, dd_div, dd_sqrt
   public :: operator(+), operator(-), operator(*), operator(/)
   public :: operator(<), operator(>), operator(<=), operator(>=), operator(==)
   public :: sqrt, abs
   public :: assignment(=), operator(**)

   interface assignment(=)
      module procedure dd_assign_r
   end interface
   interface operator(**)
      module procedure dd_pow_i
   end interface
   public :: dd_exp, dd_erf, dd_erfc
   public :: exp, erf, erfc

   interface exp
      module procedure dd_exp
   end interface
   interface erf
      module procedure dd_erf
   end interface
   interface erfc
      module procedure dd_erfc
   end interface

   ! Split so `hi + lo` reproduces the constant to ~1e-33. Hard-coded rather
   ! than derived from real128 at startup, because real128 is precisely what
   ! this module exists to avoid needing.
   !
   ! **`lo` is the remainder against the exact binary value of `hi`, not against
   ! its printed form.** Generating these from `repr(hi)` -- the shortest decimal
   ! that round-trips -- leaves `lo` wrong by the gap between the two, about
   ! 2e-17. The first version did exactly that, and the symptom was oddly
   ! specific: exp was accurate to 1e-29 while the range reduction chose k = 0
   ! and fell to 1.4e-17 the moment it chose anything else, because that is when
   ! LN2 first enters the arithmetic.
   type(dd), parameter :: DD_LN2 = dd(0.6931471805599453_dp, 2.3190468138462996e-17_dp)
   type(dd), parameter :: DD_INV_LN2 = dd(1.4426950408889634_dp, 2.0355273740931033e-17_dp)
   type(dd), parameter :: DD_2_SQRTPI = dd(1.1283791670955126_dp, 1.533545961316588e-17_dp)

   ! ERF_MAC: 23 terms, every pair verified to <1e-31
   type(dd), parameter :: ERF_MAC(0:22) = [ &
      dd(1.1283791670955126_dp, 1.533545961316588e-17_dp), &
      dd(-0.37612638903183754_dp, 1.3391897206030649e-17_dp), &
      dd(0.11283791670955126_dp, -4.017569161809194e-18_dp), &
      dd(-0.026866170645131252_dp, 4.6092880729453e-19_dp), &
      dd(0.005223977625442188_dp, -8.962504586282528e-20_dp), &
      dd(-0.0008548327023450853_dp, 5.0148896786169737e-20_dp), &
      dd(0.00012055332981789664_dp, 6.480246840070509e-21_dp), &
      dd(-1.492565035840625e-05_dp, -6.248427055364001e-22_dp), &
      dd(1.6462114365889248e-06_dp, -1.0547266132407653e-22_dp), &
      dd(-1.6365844691234924e-07_dp, 1.5075323135139275e-24_dp), &
      dd(1.4807192815879218e-08_dp, -3.254656350443331e-25_dp), &
      dd(-1.2290555301717928e-09_dp, 9.976105519856072e-26_dp), &
      dd(9.422759064650411e-11_dp, -2.8231273253303265e-27_dp), &
      dd(-6.7113668551641105e-12_dp, 1.0441838553137576e-28_dp), &
      dd(4.4632242632864775e-13_dp, -1.415652238799887e-29_dp), &
      dd(-2.7835162072109215e-14_dp, 1.4189660114017355e-30_dp), &
      dd(1.6342614095367152e-15_dp, -1.159587670993415e-32_dp), &
      dd(-9.063970842808673e-17_dp, 3.004743561097066e-33_dp), &
      dd(4.763348040515068e-18_dp, -1.856680962530252e-34_dp), &
      dd(-2.3784598852774293e-19_dp, -1.5875374697145176e-35_dp), &
      dd(1.131218725924631e-20_dp, 7.245880865990418e-38_dp), &
      dd(-5.136209054585811e-22_dp, -1.4537189752115876e-38_dp), &
      dd(2.2308786802746453e-23_dp, -4.7283897669065125e-40_dp) ]

   ! CXB: 35 terms, every pair verified to <1e-31
   type(dd), parameter :: CXB(0:34) = [ &
      dd(0.4062917348653128_dp, 1.6774286231500304e-17_dp), &
      dd(-0.18151089830658812_dp, 2.8324875440536554e-18_dp), &
      dd(0.036282977131463424_dp, -5.19682203063749e-19_dp), &
      dd(-0.0066414453256516195_dp, -4.137921097469651e-19_dp), &
      dd(0.0011300743535714113_dp, 2.4209042235881436e-20_dp), &
      dd(-0.0001806468989853791_dp, 5.023431378875334e-21_dp), &
      dd(2.7342880285563356e-05_dp, -1.1536720979841184e-21_dp), &
      dd(-3.942598196325111e-06_dp, -1.7419782678376548e-22_dp), &
      dd(5.441792670106758e-07_dp, 4.3459190761630427e-23_dp), &
      dd(-7.218251070104982e-08_dp, -3.991412939068519e-24_dp), &
      dd(9.231522167631464e-09_dp, -4.580877483699805e-25_dp), &
      dd(-1.1414743168430211e-09_dp, 6.461005232597751e-26_dp), &
      dd(1.367848537001346e-10_dp, 8.118746762064914e-27_dp), &
      dd(-1.591764656641901e-11_dp, 2.9914364538568396e-28_dp), &
      dd(1.802051255167746e-12_dp, -7.603268739698745e-29_dp), &
      dd(-1.9878778583973728e-13_dp, -9.170148516580802e-30_dp), &
      dd(2.1397210370699588e-14_dp, -1.1486646817682185e-30_dp), &
      dd(-2.2501737469358202e-15_dp, -1.2232136570604997e-31_dp), &
      dd(2.3145069843923556e-16_dp, 1.897108841112813e-32_dp), &
      dd(-2.3309339983129675e-17_dp, 1.874019689098029e-34_dp), &
      dd(2.300571581031184e-18_dp, -7.201050728147661e-35_dp), &
      dd(-2.227132705200476e-19_dp, 3.1205341724755424e-38_dp), &
      dd(2.1164143036953425e-20_dp, 1.3860885668884489e-36_dp), &
      dd(-1.9756603241087112e-21_dp, 1.4944673632591166e-37_dp), &
      dd(1.8128830779636716e-22_dp, -4.825825454946206e-39_dp), &
      dd(-1.6362173358385416e-23_dp, -5.996760949246596e-40_dp), &
      dd(1.4533650095800887e-24_dp, 6.063035670947172e-42_dp), &
      dd(-1.2711675650832093e-25_dp, -5.078494822192844e-42_dp), &
      dd(1.0953222350192672e-26_dp, -3.5199034868300316e-43_dp), &
      dd(-9.30239572291492e-28_dp, -6.364304823177461e-44_dp), &
      dd(7.790257611782451e-29_dp, -1.0256863555859985e-45_dp), &
      dd(-6.435641585129341e-30_dp, 6.287428151470075e-46_dp), &
      dd(5.246666253200263e-31_dp, -3.245655691359941e-47_dp), &
      dd(-4.222655315980911e-32_dp, 1.0771661635863134e-48_dp), &
      dd(3.356207203483348e-33_dp, 3.1244013848598736e-49_dp) ]

   ! CXC: 55 terms, every pair verified to <1e-31
   type(dd), parameter :: CXC(0:54) = [ &
      dd(0.12688864224726762_dp, 6.096112571785444e-18_dp), &
      dd(-0.08616200556789433_dp, 1.9525289463403056e-18_dp), &
      dd(0.028651053925305044_dp, -1.0724557642268018e-18_dp), &
      dd(-0.009343253046868725_dp, 8.210245882909174e-19_dp), &
      dd(0.0029914732487633434_dp, -8.923089771810032e-20_dp), &
      dd(-0.0009413219492765497_dp, -3.9559059104691374e-20_dp), &
      dd(0.0002913695785780654_dp, -1.0677164842276307e-21_dp), &
      dd(-8.878685953005519e-05_dp, -1.2701922743970771e-21_dp), &
      dd(2.6653991769316863e-05_dp, -1.6379342055070324e-21_dp), &
      dd(-7.887979843277031e-06_dp, 7.121423128999653e-22_dp), &
      dd(2.3025797005119897e-06_dp, -6.772738881748809e-23_dp), &
      dd(-6.633491767511841e-07_dp, -2.5831546778501326e-23_dp), &
      dd(1.886958777670757e-07_dp, 6.5373445953276504e-24_dp), &
      dd(-5.3024001301092376e-08_dp, -1.699429663716183e-24_dp), &
      dd(1.472498041089986e-08_dp, 6.839847001650318e-25_dp), &
      dd(-4.042761890680327e-09_dp, 1.8941433719554816e-25_dp), &
      dd(1.0977415380095793e-09_dp, 1.2475571165348304e-26_dp), &
      dd(-2.948949070772245e-10_dp, -1.802669604524215e-26_dp), &
      dd(7.840015906100578e-11_dp, 6.399200391900076e-27_dp), &
      dd(-2.0633733117411154e-11_dp, 1.3696865200671105e-27_dp), &
      dd(5.377380977005958e-12_dp, 3.681924529764916e-29_dp), &
      dd(-1.3880650612558013e-12_dp, -6.215359053087517e-29_dp), &
      dd(3.5497871470216253e-13_dp, -4.673795806694174e-31_dp), &
      dd(-8.99600830079053e-14_dp, -3.622539201200953e-31_dp), &
      dd(2.259688643949653e-14_dp, 6.5210406597575405e-31_dp), &
      dd(-5.627166849574083e-15_dp, -9.194990545971325e-32_dp), &
      dd(1.3895048089718878e-15_dp, 7.623956382873621e-32_dp), &
      dd(-3.402842116683346e-16_dp, 1.4770039164568218e-32_dp), &
      dd(8.266344962636707e-17_dp, -4.415142842117253e-33_dp), &
      dd(-1.992278713386956e-17_dp, -1.2996166017484405e-33_dp), &
      dd(4.764562143034114e-18_dp, 2.268711800854185e-34_dp), &
      dd(-1.1308384542322735e-18_dp, -2.536028064459608e-35_dp), &
      dd(2.6640858570727695e-19_dp, 1.0059960525788455e-35_dp), &
      dd(-6.230581799329958e-20_dp, 1.944129704378642e-36_dp), &
      dd(1.4467777281769567e-20_dp, -6.11950203038403e-37_dp), &
      dd(-3.335996114738264e-21_dp, 5.458501083031246e-38_dp), &
      dd(7.639328513634148e-22_dp, -9.295238706310079e-39_dp), &
      dd(-1.7375769923301635e-22_dp, 3.482765204452262e-39_dp), &
      dd(3.9259380127289167e-23_dp, -1.6048481611351936e-40_dp), &
      dd(-8.812573037525616e-24_dp, 3.389742917957844e-40_dp), &
      dd(1.9654841633846056e-24_dp, 1.4081903068057502e-41_dp), &
      dd(-4.356019646633917e-25_dp, 3.7327668115650966e-41_dp), &
      dd(9.594181088393498e-26_dp, -2.3189348250990255e-42_dp), &
      dd(-2.1002288654274108e-26_dp, 1.4697630641453873e-43_dp), &
      dd(4.569907008967642e-27_dp, -3.8523599497099677e-44_dp), &
      dd(-9.884851977740967e-28_dp, -8.223967100177135e-44_dp), &
      dd(2.1256599704494173e-28_dp, 5.673899828989284e-46_dp), &
      dd(-4.544810254409559e-29_dp, -2.0839826268359592e-45_dp), &
      dd(9.662119621544905e-30_dp, 1.4194433662743134e-46_dp), &
      dd(-2.042673696891575e-30_dp, 1.4195722915996438e-46_dp), &
      dd(4.294668991731172e-31_dp, 4.301198376361296e-47_dp), &
      dd(-8.980441293737488e-32_dp, 1.7993125474452967e-48_dp), &
      dd(1.8678205255238366e-32_dp, 5.8415455231002554e-49_dp), &
      dd(-3.864322569078916e-33_dp, 3.2732264522113656e-49_dp), &
      dd(7.95321562979717e-34_dp, 1.5994817461838743e-50_dp) ]

   ! ASY_D: 56 terms; largest |coef| = 9.64e+71
   type(dd), parameter :: ASY_D(0:55) = [ &
      dd(1.0_dp, 0.0_dp), &
      dd(-0.5_dp, 0.0_dp), &
      dd(0.75_dp, 0.0_dp), &
      dd(-1.875_dp, 0.0_dp), &
      dd(6.5625_dp, 0.0_dp), &
      dd(-29.53125_dp, 0.0_dp), &
      dd(162.421875_dp, 0.0_dp), &
      dd(-1055.7421875_dp, 0.0_dp), &
      dd(7918.06640625_dp, 0.0_dp), &
      dd(-67303.564453125_dp, 0.0_dp), &
      dd(639383.8623046875_dp, 0.0_dp), &
      dd(-6713530.554199219_dp, 0.0_dp), &
      dd(77205601.37329102_dp, 0.0_dp), &
      dd(-965070017.1661377_dp, 0.0_dp), &
      dd(13028445231.742859_dp, 0.0_dp), &
      dd(-188912455860.27145_dp, 0.0_dp), &
      dd(2928143065834.2075_dp, 1.52587890625e-05_dp), &
      dd(-48314360586264.42_dp, -0.00244903564453125_dp), &
      dd(845501310259627.4_dp, 0.050670623779296875_dp), &
      dd(-1.5641774239803108e+16_dp, 0.6250934600830078_dp), &
      dd(3.050145976761606e+17_dp, 17.810677528381348_dp), &
      dd(-6.252799252361292e+18_dp, -397.1188893318176_dp), &
      dd(1.3443518392576778e+20_dp, -677.943879365921_dp), &
      dd(-3.024791638329775e+21_dp, -255082.26271426678_dp), &
      dd(7.108260350074972e+22_dp, -2918462.8262147307_dp), &
      dd(-1.741523785768368e+24_dp, -37549564.7577391_dp), &
      dd(4.440885653709338e+25_dp, 2299691181.322347_dp), &
      dd(-1.1768346982329746e+27_dp, -60941816305.0422_dp), &
      dd(3.2362954201406804e+28_dp, -110806446747.33965_dp), &
      dd(-9.223441947400939e+29_dp, -45220527889844.82_dp), &
      dd(2.720915374483277e+31_dp, -1832587915244457.8_dp), &
      dd(-8.298791892173995e+32_dp, 5.364213160127071e+16_dp), &
      dd(2.6141194460348083e+34_dp, 1.2646342101150177e+18_dp), &
      dd(-8.495888199613127e+35_dp, 2.710405346322105e+18_dp), &
      dd(2.8461225468703976e+37_dp, 8.684321127311062e+20_dp), &
      dd(-9.819122786702872e+38_dp, 1.0179207215168822e+22_dp), &
      dd(3.4857885892795196e+40_dp, -5.8803544731623615e+23_dp), &
      dd(-1.2723128350870246e+42_dp, -2.6893738957542547e+25_dp), &
      dd(4.771173131576342e+43_dp, 3.40702403702327e+27_dp), &
      dd(-1.8369016556568918e+45_dp, 1.7382379288849744e+28_dp), &
      dd(7.255761539844723e+46_dp, -5.27786817954502e+28_dp), &
      dd(-2.9385834236371126e+48_dp, -8.003668189110103e+30_dp), &
      dd(1.2195121208094018e+50_dp, -1.6149590921024911e+33_dp), &
      dd(-5.1829265134399576e+51_dp, 1.9325088601919174e+35_dp), &
      dd(2.2545730333463817e+53_dp, -1.9040237508114167e+37_dp), &
      dd(-1.0032849998391398e+55_dp, 4.857405542575833e+38_dp), &
      dd(4.564946749268086e+56_dp, 1.6690994610266944e+40_dp), &
      dd(-2.12270023840966e+58_dp, 8.354460403601517e+41_dp), &
      dd(1.0082826132445884e+60_dp, 3.279373497811732e+43_dp), &
      dd(-4.890170674236254e+61_dp, 1.5316081813555972e+45_dp), &
      dd(2.4206344837469458e+63_dp, 3.265621966855089e+46_dp), &
      dd(-1.2224204142922075e+65_dp, -9.504710393915423e+48_dp), &
      dd(6.295465133604869e+66_dp, 4.519608753804979e+49_dp), &
      dd(-3.305119195142556e+68_dp, -1.4345416008762372e+52_dp), &
      dd(1.7682387694012676e+70_dp, 7.914249992948164e+53_dp), &
      dd(-9.636901293236909e+71_dp, -1.8612733807713268e+55_dp) ]

   ! **Overloaded so the shared bodies do not change.** cint_fmt_body.inc and
   ! cint_wheeler_body.inc are included once per precision with `rk` bound to a
   ! kind; binding a derived type instead only works if `a*b - c/d` still reads
   ! as arithmetic. Everything below exists so that the algorithm text stays
   ! one copy rather than becoming two that must be kept in step.
   interface operator(+)
      module procedure dd_add, dd_add_r, dd_r_add
   end interface
   interface operator(-)
      module procedure dd_sub, dd_sub_r, dd_r_sub, dd_neg
   end interface
   interface operator(*)
      module procedure dd_mul, dd_mul_r, dd_r_mul
   end interface
   interface operator(/)
      module procedure dd_div, dd_div_r, dd_r_div
   end interface
   interface operator(<)
      module procedure dd_lt, dd_lt_r
   end interface
   interface operator(>)
      module procedure dd_gt, dd_gt_r
   end interface
   interface operator(<=)
      module procedure dd_le, dd_le_r
   end interface
   interface operator(>=)
      module procedure dd_ge, dd_ge_r
   end interface
   interface operator(==)
      module procedure dd_eq, dd_eq_r
   end interface
   interface sqrt
      module procedure dd_sqrt
   end interface
   interface abs
      module procedure dd_abs
   end interface

   ! Dekker's splitting constant, 2**27 + 1: splits a 53-bit significand into
   ! two 26-bit halves whose product is exact.
   real(dp), parameter :: SPLITTER = 134217729.0_dp

contains

   pure type(dd) function dd_from(x) result(r)
      !! A double promoted exactly: the error term is zero by construction.
      real(dp), intent(in) :: x
      r%hi = x
      r%lo = 0.0_dp
   end function dd_from

   pure real(dp) function dd_to_dp(a) result(r)
      !! Back to double, correctly rounded by the representation itself.
      type(dd), intent(in) :: a
      r = a%hi + a%lo
   end function dd_to_dp

   pure subroutine two_sum(a, b, s, e)
      !! s = fl(a+b) and e = the exact discarded error, so s + e == a + b.
      !!
      !! Knuth's version: no assumption about which operand is larger, at the
      !! cost of three more flops than Dekker's fast_two_sum.
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: s, e
      real(dp) :: bb

      s = a + b
      bb = s - a
      e = (a - (s - bb)) + (b - bb)
   end subroutine two_sum

   pure subroutine split(a, hi, lo)
      !! a = hi + lo exactly, each with a 26-bit significand.
      real(dp), intent(in) :: a
      real(dp), intent(out) :: hi, lo
      real(dp) :: t

      t = SPLITTER * a
      hi = t - (t - a)
      lo = a - hi
   end subroutine split

   pure subroutine two_prod(a, b, p, e)
      !! p = fl(a*b) and e = the exact discarded error, so p + e == a * b.
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: p, e
      real(dp) :: ahi, alo, bhi, blo

      p = a * b
      call split(a, ahi, alo)
      call split(b, bhi, blo)
      e = ((ahi * bhi - p) + ahi * blo + alo * bhi) + alo * blo
   end subroutine two_prod

   pure type(dd) function dd_add(a, b) result(r)
      !! Sum of two double-doubles, renormalised.
      type(dd), intent(in) :: a, b
      real(dp) :: s1, s2, t1, t2, u1, u2

      ! Every `two_sum` writes to variables it does not also read. Passing the
      ! same name as both an intent(in) and an intent(out) argument is aliasing,
      ! which Fortran leaves undefined -- and the first version of this routine
      ! did exactly that, renormalising through `two_sum(s1, s2, s1, s2)`. It
      ! compiled without a murmur and returned a sum of roughly zero.
      call two_sum(a%hi, b%hi, s1, s2)
      call two_sum(a%lo, b%lo, t1, t2)
      s2 = s2 + t1
      call two_sum(s1, s2, u1, u2)
      u2 = u2 + t2
      call two_sum(u1, u2, r%hi, r%lo)
   end function dd_add

   pure type(dd) function dd_mul(a, b) result(r)
      !! Product of two double-doubles, renormalised.
      type(dd), intent(in) :: a, b
      real(dp) :: p1, p2

      call two_prod(a%hi, b%hi, p1, p2)
      p2 = p2 + a%hi * b%lo + a%lo * b%hi
      call two_sum(p1, p2, r%hi, r%lo)
   end function dd_mul

   pure type(dd) function dd_neg(a) result(r)
      type(dd), intent(in) :: a
      r%hi = -a%hi; r%lo = -a%lo
   end function dd_neg

   pure type(dd) function dd_sub(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = dd_add(a, dd_neg(b))
   end function dd_sub

   pure type(dd) function dd_div(a, b) result(r)
      !! Long division: a double-precision quotient, then one correction pass
      !! against the exact remainder.
      !!
      !! The remainder is what makes it work. `a - q1*b` is computed in
      !! double-double, so the leading terms cancel exactly rather than to
      !! within a rounding, and what survives is the part the first quotient
      !! missed.
      type(dd), intent(in) :: a, b
      real(dp) :: q1, q2
      type(dd) :: r1, prod

      q1 = a%hi / b%hi
      prod = dd_mul(b, dd_from(q1))
      r1 = dd_sub(a, prod)
      q2 = (r1%hi + r1%lo) / b%hi
      call two_sum(q1, q2, r%hi, r%lo)
   end function dd_div

   pure type(dd) function dd_sqrt(a) result(r)
      !! One Newton step on a double-precision square root.
      !!
      !! x1 = x0 + (a - x0^2) / (2 x0) doubles the correct digits, so a
      !! 53-bit seed lands at the ~106 the type carries. The residual `a - x0^2`
      !! is formed in double-double for the same reason as in the division:
      !! computed in double it would be pure rounding noise.
      type(dd), intent(in) :: a
      real(dp) :: x0
      type(dd) :: resid, corr

      if (a%hi <= 0.0_dp) then
         r%hi = 0.0_dp; r%lo = 0.0_dp
         return
      end if
      x0 = sqrt(a%hi)
      resid = dd_sub(a, dd_mul(dd_from(x0), dd_from(x0)))
      corr = dd_div(resid, dd_from(2.0_dp * x0))
      r = dd_add(dd_from(x0), corr)
   end function dd_sqrt

   pure type(dd) function dd_abs(a) result(r)
      type(dd), intent(in) :: a
      if (a%hi < 0.0_dp) then
         r = dd_neg(a)
      else
         r = a
      end if
   end function dd_abs

   ! Mixed dd/real64 forms, so a literal in the shared bodies needs no wrapper.
   pure type(dd) function dd_add_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_add(a, dd_from(b))
   end function dd_add_r

   pure type(dd) function dd_r_add(a, b) result(r)
      real(dp), intent(in) :: a
      type(dd), intent(in) :: b
      r = dd_add(dd_from(a), b)
   end function dd_r_add

   pure type(dd) function dd_sub_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_sub(a, dd_from(b))
   end function dd_sub_r

   pure type(dd) function dd_r_sub(a, b) result(r)
      real(dp), intent(in) :: a
      type(dd), intent(in) :: b
      r = dd_sub(dd_from(a), b)
   end function dd_r_sub

   pure type(dd) function dd_mul_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_mul(a, dd_from(b))
   end function dd_mul_r

   pure type(dd) function dd_r_mul(a, b) result(r)
      real(dp), intent(in) :: a
      type(dd), intent(in) :: b
      r = dd_mul(dd_from(a), b)
   end function dd_r_mul

   pure type(dd) function dd_div_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_div(a, dd_from(b))
   end function dd_div_r

   pure type(dd) function dd_r_div(a, b) result(r)
      real(dp), intent(in) :: a
      type(dd), intent(in) :: b
      r = dd_div(dd_from(a), b)
   end function dd_r_div

   ! Comparisons go through the full value, not just `hi`: two numbers can share
   ! a leading double and differ in the tail, which is the entire reason for the
   ! type.
   pure logical function dd_lt(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = (a%hi < b%hi) .or. (a%hi == b%hi .and. a%lo < b%lo)
   end function dd_lt

   pure logical function dd_lt_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_lt(a, dd_from(b))
   end function dd_lt_r

   pure logical function dd_gt(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = (a%hi > b%hi) .or. (a%hi == b%hi .and. a%lo > b%lo)
   end function dd_gt

   pure logical function dd_gt_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_gt(a, dd_from(b))
   end function dd_gt_r

   pure logical function dd_le(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = .not. dd_gt(a, b)
   end function dd_le

   pure logical function dd_le_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_le(a, dd_from(b))
   end function dd_le_r

   pure logical function dd_ge(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = .not. dd_lt(a, b)
   end function dd_ge

   pure logical function dd_ge_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_ge(a, dd_from(b))
   end function dd_ge_r

   pure logical function dd_eq(a, b) result(r)
      type(dd), intent(in) :: a, b
      r = (a%hi == b%hi) .and. (a%lo == b%lo)
   end function dd_eq

   pure logical function dd_eq_r(a, b) result(r)
      type(dd), intent(in) :: a
      real(dp), intent(in) :: b
      r = dd_eq(a, dd_from(b))
   end function dd_eq_r

   pure type(dd) function dd_exp(a) result(r)
      !! exp for double-double, by range reduction and a Taylor series.
      !!
      !! Two reductions. First `a = k ln2 + t` puts `t` in [-ln2/2, ln2/2]; then
      !! `t` is halved NSQ more times and the result squared back up.
      !!
      !! **NSQ is a trade, and it was measured rather than guessed.** Each
      !! squaring doubles the relative error, so NSQ of 16 puts a floor of
      !! 2^16 * 1e-32 under the answer -- 2.3e-27, measured. Fewer squarings
      !! mean a larger series argument and more terms. Over [-35, 35]: NSQ 16
      !! gives 2.3e-27, NSQ 8 gives 1.4e-29, NSQ 6 gives 3.3e-30. Below 6 the
      !! twelve terms below stop converging and it gets worse again.
      type(dd), intent(in) :: a
      integer, parameter :: NSQ = 6
      type(dd) :: t, term, sum, sq
      real(dp) :: k
      integer :: i

      ! ln(tiny(1.0_dp)) = -708.396, so anything above this still has a normal
      ! double for `hi`. The first version cut at a round -700, which is
      ! exp = 9.8e-305 -- nowhere near underflow, and it silently returned zero
      ! for a whole decade of representable results. erfc(26.46) came back 0
      ! against a true 1.8e-306.
      if (a%hi < -708.0_dp) then
         r = dd_from(0.0_dp)
         return
      end if

      k = anint(dd_to_dp(dd_mul(a, DD_INV_LN2)))
      t = dd_sub(a, dd_mul(DD_LN2, dd_from(k)))
      ! Halve NSQ times: exp(t) = exp(t/2^NSQ)^(2^NSQ)
      t = dd_mul(t, dd_from(2.0_dp**(-NSQ)))

      ! 1 + t + t^2/2! + ...
      sum = dd_add(dd_from(1.0_dp), t)
      term = t
      do i = 2, 12
         term = dd_mul(term, t)
         term = dd_div(term, dd_from(real(i, dp)))
         sum = dd_add(sum, term)
      end do

      do i = 1, NSQ
         sum = dd_mul(sum, sum)
      end do

      r = dd_mul(sum, dd_from(2.0_dp**int(k)))
   end function dd_exp

   pure type(dd) function dd_erf(a) result(r)
      !! erf for double-double.
      !!
      !! Below the crossover this uses the *confluent* series
      !!
      !!     erf(x) = (2x/sqrt(pi)) exp(-x^2) sum_n (2x^2)^n / (1.3.5...(2n+1))
      !!
      !! rather than the textbook alternating one. Every term is positive, so
      !! there is no cancellation; the alternating series loses roughly x^2/ln10
      !! digits to it and is useless here by x = 4, which is inside the range
      !! the Boys function asks for.
      type(dd), intent(in) :: a
      type(dd) :: x, x2, term, sum, two_x2
      real(dp) :: den
      integer :: n

      x = dd_abs(a)
      if (x%hi > 6.5_dp) then
         r = dd_from(1.0_dp)
         if (a%hi < 0.0_dp) r = dd_neg(r)
         return
      end if

      x2 = dd_mul(x, x)
      two_x2 = dd_mul(x2, dd_from(2.0_dp))
      term = dd_from(1.0_dp)
      sum = dd_from(1.0_dp)
      do n = 1, 200
         den = real(2*n + 1, dp)
         term = dd_div(dd_mul(term, two_x2), dd_from(den))
         sum = dd_add(sum, term)
         if (abs(term%hi) < 1.0e-36_dp * abs(sum%hi)) exit
      end do

      r = dd_mul(dd_mul(DD_2_SQRTPI, x), dd_mul(dd_exp(dd_neg(x2)), sum))
      if (a%hi < 0.0_dp) r = dd_neg(r)
   end function dd_erf

   pure type(dd) function dd_erfc(a) result(r)
      !! erfc for double-double, in four bands.
      !!
      !! Band structure and method from bitwise_adventures' `erfc_reprod`, which
      !! is the same shape at double precision. Two things had to change for dd.
      !!
      !! **The C/D boundary moves from 7 to 9.** An asymptotic series cannot be
      !! made more accurate than its smallest term, and for this one that floor
      !! is 7.4e-22 at x = 7 -- perfect for double and unreachable for dd. At
      !! x = 9 it is 9.6e-36. So the Chebyshev band has to cover [2, 9] and the
      !! asymptotic band start where it can actually deliver.
      !!
      !! **`1 - erf(x)` is used only below 0.46875**, where erf is small enough
      !! that the subtraction costs nothing. Above it erfc is computed directly,
      !! because by x = 4 the identity throws away eight digits before any of
      !! this arithmetic begins -- and the two call sites in the shared bodies
      !! are `erfc(a) - erfc(b)` with the two nearly equal, which is where those
      !! digits were going to be needed.
      type(dd), intent(in) :: a
      type(dd) :: ax, u, t, tt, b1, b2, bt, ex2, s2
      integer :: j

      ax = dd_abs(a)
      if (ax%hi >= 28.0_dp) then
         ! Underflowed; erfc(28) is ~1e-342, below what a double `hi` can hold.
         r = dd_from(0.0_dp)
         if (a%hi < 0.0_dp) r = dd_from(2.0_dp)
         return
      end if

      if (ax%hi <= 0.46875_dp) then
         u = dd_mul(a, a)
         t = ERF_MAC(ubound(ERF_MAC, 1))
         do j = ubound(ERF_MAC, 1) - 1, 0, -1
            t = dd_add(ERF_MAC(j), dd_mul(u, t))
         end do
         r = dd_sub(dd_from(1.0_dp), dd_mul(a, t))
         return
      end if

      ! erfcx(x) = erfc(x) exp(x^2), which is what the fits below approximate;
      ! it is O(1) across the whole range where erfc itself spans 300 decades.
      ex2 = dd_exp(dd_neg(dd_mul(ax, ax)))

      if (ax%hi <= 2.0_dp) then
         ! Clenshaw on [0.46875, 2]. Both constants are exact in binary.
         tt = dd_div(dd_sub(dd_mul(ax, dd_from(2.0_dp)), dd_from(2.46875_dp)), &
                     dd_from(1.53125_dp))
         b1 = dd_from(0.0_dp); b2 = dd_from(0.0_dp)
         do j = ubound(CXB, 1), 1, -1
            bt = dd_add(dd_sub(dd_mul(dd_mul(dd_from(2.0_dp), tt), b1), b2), CXB(j))
            b2 = b1; b1 = bt
         end do
         t = dd_add(dd_sub(dd_mul(tt, b1), b2), CXB(0))
         r = dd_mul(ex2, t)
      else if (ax%hi <= 9.0_dp) then
         tt = dd_div(dd_sub(dd_mul(ax, dd_from(2.0_dp)), dd_from(11.0_dp)), &
                     dd_from(7.0_dp))
         b1 = dd_from(0.0_dp); b2 = dd_from(0.0_dp)
         do j = ubound(CXC, 1), 1, -1
            bt = dd_add(dd_sub(dd_mul(dd_mul(dd_from(2.0_dp), tt), b1), b2), CXC(j))
            b2 = b1; b1 = bt
         end do
         t = dd_add(dd_sub(dd_mul(tt, b1), b2), CXC(0))
         r = dd_mul(ex2, t)
      else
         ! Asymptotic in s = 1/x^2, summed backwards so the smallest terms go in
         ! first. Divergent, but cut where the terms stop shrinking.
         s2 = dd_div(dd_from(1.0_dp), dd_mul(ax, ax))
         t = ASY_D(ubound(ASY_D, 1))
         do j = ubound(ASY_D, 1) - 1, 0, -1
            t = dd_add(ASY_D(j), dd_mul(s2, t))
         end do
         r = dd_div(dd_mul(ex2, dd_mul(dd_div(DD_2_SQRTPI, dd_from(2.0_dp)), t)), ax)
      end if

      if (a%hi < 0.0_dp) r = dd_sub(dd_from(2.0_dp), r)
   end function dd_erfc

   elemental subroutine dd_assign_r(out, in)
      !! `x = 0.5_dp` where x is a dd. Needed because the shared bodies assign
      !! plain literals to working variables, and without this every one is
      !! "Cannot convert REAL(8) to TYPE(dd)".
      type(dd), intent(out) :: out
      real(dp), intent(in) :: in
      out%hi = in
      out%lo = 0.0_dp
   end subroutine dd_assign_r

   pure type(dd) function dd_pow_i(a, n) result(r)
      !! Integer power by repeated squaring.
      !!
      !! Squaring rather than exp(n log a): it is exact for the small n the
      !! bodies use, and it keeps a negative base working, which the log form
      !! does not.
      type(dd), intent(in) :: a
      integer, intent(in) :: n
      type(dd) :: base
      integer :: k

      if (n == 0) then
         r = dd_from(1.0_dp)
         return
      end if
      base = a
      k = abs(n)
      r = dd_from(1.0_dp)
      do while (k > 0)
         if (mod(k, 2) == 1) r = dd_mul(r, base)
         base = dd_mul(base, base)
         k = k/2
      end do
      if (n < 0) r = dd_div(dd_from(1.0_dp), r)
   end function dd_pow_i

end module cint_dd
