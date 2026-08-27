/* binary128 arithmetic for Fortran compilers that have no real128.
 *
 * nvfortran and LLVM Flang both report `real128` as -1, so `real(qp)` is a
 * hard error there. The precision itself is not missing: nvc and clang both
 * have __float128 with a true 113-bit significand, and libquadmath supplies
 * the transcendentals. This exposes both to Fortran through bind(c).
 *
 * **Everything is passed by reference.** A 16-byte struct by value has no
 * settled ABI across gfortran/nvfortran/flang and their C partners, and this
 * exists to be identical on all three. By reference there is no question.
 *
 * libquadmath ships with GCC, so its header and archive come from there even
 * when the Fortran driver is nvfortran or flang; both link it happily, and
 * libquadmath.a makes it static so no .so has to be found at run time.
 */
#include <quadmath.h>

typedef __float128 q_t;
#define A (*(const q_t*)a)
#define B (*(const q_t*)b)
#define R (*(q_t*)r)

void mqq_add(const void*a,const void*b,void*r){R=A+B;}
void mqq_sub(const void*a,const void*b,void*r){R=A-B;}
void mqq_mul(const void*a,const void*b,void*r){R=A*B;}
void mqq_div(const void*a,const void*b,void*r){R=A/B;}
void mqq_neg(const void*a,void*r){R=-A;}
void mqq_abs(const void*a,void*r){R=fabsq(A);}
void mqq_sqrt(const void*a,void*r){R=sqrtq(A);}
void mqq_exp(const void*a,void*r){R=expq(A);}
void mqq_erf(const void*a,void*r){R=erfq(A);}
void mqq_erfc(const void*a,void*r){R=erfcq(A);}
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
