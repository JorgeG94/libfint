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
 * The transcendentals are glibc's *f128 family, NOT libquadmath's *q family,
 * and the difference is measured, not stylistic: gfortran compiles real128
 * sqrt/exp/erf/erfc to sqrtf128/expf128/erff128/erfcf128 (nm a real128
 * object and look), and libquadmath's sqrtq disagrees with glibc's
 * correctly-rounded sqrtf128 in the last bit -- 526 of 10,000 values in
 * quad_check when this shim was first written against quadmath.h.
 * Bit-identity with gfortran's real128 is this module's whole contract, so
 * the shim calls exactly the symbols gfortran calls. They live in libm
 * (glibc 2.26+); the declarations are written out here because <math.h>
 * only offers them behind the _Float128 keyword, which not every C compiler
 * this shim targets has. __float128 and _Float128 share one ABI, so the
 * link resolves either way. The +,-,*,/ below are compiler-generated
 * (libgcc's __addtf3 family), the same code gfortran emits.
 */
typedef __float128 q_t;

extern q_t sqrtf128(q_t);
extern q_t expf128(q_t);
extern q_t erff128(q_t);
extern q_t erfcf128(q_t);
extern q_t fabsf128(q_t);
#define A (*(const q_t*)a)
#define B (*(const q_t*)b)
#define R (*(q_t*)r)

void mqq_add(const void*a,const void*b,void*r){R=A+B;}
void mqq_sub(const void*a,const void*b,void*r){R=A-B;}
void mqq_mul(const void*a,const void*b,void*r){R=A*B;}
void mqq_div(const void*a,const void*b,void*r){R=A/B;}
void mqq_neg(const void*a,void*r){R=-A;}
void mqq_abs(const void*a,void*r){R=fabsf128(A);}
void mqq_sqrt(const void*a,void*r){R=sqrtf128(A);}
void mqq_exp(const void*a,void*r){R=expf128(A);}
void mqq_erf(const void*a,void*r){R=erff128(A);}
void mqq_erfc(const void*a,void*r){R=erfcf128(A);}
void mqq_from_d(double d,void*r){R=(q_t)d;}
/* Three doubles reassembled. binary128 needs 113 bits and a double holds 53,
   so two are short by seven and three are exact with room over. The sum is
   exact because the parts are non-overlapping by construction. This is how a
   67-digit constant reaches a type Fortran cannot write a literal for. */
void mqq_from_3d(double a1,double a2,double a3,void*r){R=((q_t)a1+(q_t)a2)+(q_t)a3;}
void mqq_from_i(int n,void*r){R=(q_t)n;}
double mqq_to_d(const void*a){return (double)A;}
/* Comparison returns int rather than a Fortran logical: the bit pattern of
   .true. is not fixed across compilers, but 0/1 is. */
int mqq_cmp(const void*a,const void*b){return A<B?-1:(A>B?1:0);}
