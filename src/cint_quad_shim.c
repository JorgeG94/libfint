/* binary128 arithmetic for Fortran compilers that have no real128.
 *
 * nvfortran and LLVM Flang both report `real128` as -1, so `real(qp)` is a
 * hard error there. The precision itself is not missing: nvc and clang both
 * have __float128 with a true 113-bit significand. This exposes it to
 * Fortran through bind(c).
 *
 * **Everything is passed by reference.** A 16-byte struct by value has no
 * settled ABI across gfortran/nvfortran/flang and their C partners, and this
 * exists to be identical on all three. By reference there is no question.
 *
 * Bit-identity with gfortran's real128 is this module's whole contract, so
 * the transcendentals are whatever symbols gfortran's real128 intrinsics
 * call ON THIS PLATFORM, and the difference is measured, not stylistic. On
 * Linux/glibc that is the *f128 family (nm a real128 object and look), and
 * libquadmath's sqrtq disagrees with glibc's correctly-rounded sqrtf128 in
 * the last bit -- 526 of 10,000 values in quad_check when this shim was
 * first written against quadmath.h. Where the C library has no *f128
 * (macOS), gfortran calls libquadmath, and so does the shim. The +,-,*,/
 * below are compiler-generated (libgcc's __addtf3 family), the same code
 * gfortran emits.
 *
 * But __float128 is a per-target property, not a per-compiler one, and the
 * assumption that every C compiler on a machine has it does not survive
 * arm64 macOS: Apple clang there rejects the keyword outright ("__float128
 * is not supported on this target"), ships no <quadmath.h>, and has a
 * `long double` that is plain binary64. The same machine's Homebrew GCC has
 * all three -- and it is the compiler that supplies gfortran, whose real128
 * this file is defined against. So the choice of C compiler is not free on
 * macOS: it has to be the GCC one.
 *
 * Hence the ladder below. It takes the widest IEEE type the compiler
 * actually offers, and where that is not 113 bits it stops the build rather
 * than hand the extended Rys ladder a silently narrower one -- cint_fmt_qp,
 * cint_schmidt and cint_wheeler are compiled against this unconditionally,
 * so a quiet fallback would not be a missing feature, it would be wrong
 * numbers in integrals that still compute.
 */

/* ------------------------------------------------------------ backend ---- */
#if defined(__SIZEOF_FLOAT128__) \
    || (defined(__FLT128_MANT_DIG__) && __FLT128_MANT_DIG__ == 113)
/* __float128 exists. Both spellings are tested because the compilers do not
   agree on which to advertise: gcc on x86 defines __SIZEOF_FLOAT128__, gcc
   on aarch64-darwin defines only the __FLT128_* set, and clang defines both
   wherever it supports the type at all. */
#define CINT_QUAD_FLOAT128 1
#elif defined(__LDBL_MANT_DIG__) && __LDBL_MANT_DIG__ == 113
/* No __float128 keyword, but `long double` IS binary128 -- aarch64 and
   s390x under clang, where the two would name one type anyway. Same 113
   bits, same correctly-rounded +,-,*,/ and sqrt, and on glibc sqrtl and
   friends ARE the *f128 entry points, so the contract holds unchanged. */
#define CINT_QUAD_LONG_DOUBLE 1
#elif defined(CINT_QUAD_ALLOW_NARROW)
/* Explicitly asked for, and deliberately loud. */
#define CINT_QUAD_LONG_DOUBLE 1
#define CINT_QUAD_NARROW 1
#warning "CINT_QUAD_ALLOW_NARROW: no binary128 on this target, so cint_quad falls back to long double -- 53 bits on arm64 macOS, 64 on x86, not 113. The extended Rys ladder is then no wider than the double one and quad_check will fail. Only for builds that exercise the double path alone."
#else
#error "cint_quad_shim.c: this C compiler has no binary128 type -- no __float128 on this target, and long double is not IEEE quad. cint_quad has no other source of 113 bits, and cint_fmt_qp/cint_schmidt/cint_wheeler are built against it unconditionally. Build the C sources with a compiler that has one -- on macOS that means the GCC that provides gfortran, not Apple clang: cmake -DCMAKE_C_COMPILER=gcc-15, or fpm --c-compiler gcc-15. Or define CINT_QUAD_ALLOW_NARROW to fall back to long double and give up the extended ladder's precision."
#endif

/* <quadmath.h> ships with GCC, not with the target's C library, so its
   presence tracks the compiler rather than the platform. */
#if defined(__has_include)
#if __has_include(<quadmath.h>)
#define CINT_QUAD_HAVE_QUADMATH_H 1
#endif
#elif defined(__GNUC__) && !defined(__clang__)
#define CINT_QUAD_HAVE_QUADMATH_H 1
#endif

#ifdef CINT_QUAD_FLOAT128
typedef __float128 q_t;

#ifdef __linux__
/* glibc's *f128 family; declared by hand because <math.h> only offers them
   behind the _Float128 keyword, which not every targeted C compiler has.
   __float128 and _Float128 share one ABI, so the link resolves either way. */
extern q_t sqrtf128(q_t);
extern q_t expf128(q_t);
extern q_t erff128(q_t);
extern q_t erfcf128(q_t);
extern q_t fabsf128(q_t);
#elif defined(CINT_QUAD_HAVE_QUADMATH_H)
/* No *f128 in the C library (macOS): gfortran's real128 intrinsics call
   libquadmath there, so matching gfortran means calling libquadmath too. */
#include <quadmath.h>
#define sqrtf128 sqrtq
#define expf128  expq
#define erff128  erfq
#define erfcf128 erfcq
#define fabsf128 fabsq
#else
#error "cint_quad_shim.c: this C compiler has __float128 but no library to evaluate it with -- neither glibc's *f128 family nor GCC's <quadmath.h>. sqrt/exp/erf/erfc have nowhere to go. Build the C sources with GCC, which carries libquadmath."
#endif

#else /* CINT_QUAD_LONG_DOUBLE */
/* One IEEE type wide enough that C's own long-double math is the whole
   implementation: sqrtl/expl/erfl/erfcl, correctly rounded on +,-,*,/ and
   sqrt exactly as the *f128 ones are. */
#include <math.h>
typedef long double q_t;
#define sqrtf128 sqrtl
#define expf128  expl
#define erff128  erfl
#define erfcf128 erfcl
#define fabsf128 fabsl
#endif

/* A long double narrower than the 16 bytes the Fortran type reserves leaves
   the rest of the slot untouched, which would make the padding whatever the
   stack held.  Nothing on either side reads it -- but quad values are
   compared as raw bits in quad_const_check, and stack residue is not a
   value.  Zero the slot first; only the narrow fallback pays for it. */
#if defined(CINT_QUAD_NARROW) && defined(__SIZEOF_LONG_DOUBLE__) \
    && __SIZEOF_LONG_DOUBLE__ < 16
#define CINT_QUAD_PAD_STORE 1
#endif

/* memcpy, never a pointer cast: the Fortran side declares quad as two
   doubles, so it is 8-byte aligned, while __float128 loads and stores may
   be emitted as 16-byte-aligned SSE moves.  Casting the pointer is UB that
   gcc and nvc happened to survive; flang laid a quad on an odd 8-byte slot
   and the movaps faulted.  memcpy through a local compiles to unaligned
   moves everywhere and costs nothing. */
static inline q_t ld(const void*p){q_t v;__builtin_memcpy(&v,p,sizeof v);return v;}
#ifdef CINT_QUAD_PAD_STORE
static inline void st(void*p,q_t v){__builtin_memset(p,0,16);__builtin_memcpy(p,&v,sizeof v);}
#else
static inline void st(void*p,q_t v){__builtin_memcpy(p,&v,sizeof v);}
#endif
#define A ld(a)
#define B ld(b)

void mqq_add(const void*a,const void*b,void*r){st(r,A+B);}
void mqq_sub(const void*a,const void*b,void*r){st(r,A-B);}
void mqq_mul(const void*a,const void*b,void*r){st(r,A*B);}
void mqq_div(const void*a,const void*b,void*r){st(r,A/B);}
void mqq_neg(const void*a,void*r){st(r,-A);}
void mqq_abs(const void*a,void*r){st(r,fabsf128(A));}
void mqq_sqrt(const void*a,void*r){st(r,sqrtf128(A));}
void mqq_exp(const void*a,void*r){st(r,expf128(A));}
void mqq_erf(const void*a,void*r){st(r,erff128(A));}
void mqq_erfc(const void*a,void*r){st(r,erfcf128(A));}
void mqq_from_d(double d,void*r){st(r,(q_t)d);}
/* Three doubles reassembled. binary128 needs 113 bits and a double holds 53,
   so two are short by seven and three are exact with room over. The sum is
   exact because the parts are non-overlapping by construction. This is how a
   67-digit constant reaches a type Fortran cannot write a literal for. */
void mqq_from_3d(double a1,double a2,double a3,void*r){st(r,((q_t)a1+(q_t)a2)+(q_t)a3);}
void mqq_from_i(int n,void*r){st(r,(q_t)n);}
double mqq_to_d(const void*a){return (double)A;}
/* Comparison returns int rather than a Fortran logical: the bit pattern of
   .true. is not fixed across compilers, but 0/1 is. */
int mqq_cmp(const void*a,const void*b){return A<B?-1:(A>B?1:0);}
