/*
 * The C ABI, driven the way a C consumer drives it.
 *
 * Every other check in this suite calls the Fortran API.  That leaves the
 * layer a C program actually links -- src/cint_c_abi.f90 -- untested, and it
 * is the layer most able to fail quietly: a wrong bind(C) name is a link
 * error the suite never sees, and a wrong buffer size is a silent overwrite.
 *
 * This links libfint and drives it through libcint's own C entry points, so a
 * missing or misnamed export fails at link time rather than at a consumer's.
 *
 * The reference numbers were produced by compiling this same file against
 * upstream libcint and recording what it printed.  That is the comparison the
 * rest of the suite makes at runtime; it cannot be made at runtime here,
 * because both libraries define the same symbols and a program links one or
 * the other, never both.
 *
 * Four things are asserted, and they fail in different ways:
 *   - the dimension helpers, which a driver needs to size its own buffers;
 *   - values against libcint, to a tolerance that is really bit-identity;
 *   - the CINT2 spelling against the CINT3 one, bitwise, since they are two
 *     names for one computation;
 *   - the optimizer path against the unoptimised one, bitwise, which is what
 *     makes an optimizer a cache rather than a second implementation.
 *
 * The Hessian family -- the eleven second derivatives an analytic Hessian is
 * built from -- is driven in both angular conventions, where everything above
 * it is spherical only.  That is what the d shell in the fixture is for: below
 * l = 2 the two conventions agree value for value, so a Cartesian export could
 * be wired to its spherical twin and no basis of s and p would say so.
 */
#include <stdio.h>
#include <math.h>
#include <string.h>

int  int2e_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint2e_sph(double*, int*, int*, int, int*, int, double*, void*);
int  int1e_ovlp_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int1e_kin_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int1e_nuc_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int2e_ip1_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
void cint2e_sph_optimizer(void**, int*, int, int*, int, double*);
void CINTdel_optimizer(void**);
int  CINTcgto_spheric(int, int*);
int  CINTcgto_cart(int, int*);
int  CINTtot_cgto_spheric(int*, int);

/* The Hessian family: the second derivatives pyscf.hessian.rhf drives, in
   both angular conventions and both spellings.  Every name is written out
   rather than reached through a macro, because the name is the thing under
   test: a misspelt bind(C) is a link error only if something spells the
   symbol the way a consumer would. */
typedef int (*cint3_t)(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
typedef int (*cint2_1e_t)(double*, int*, int*, int, int*, int, double*);
typedef int (*cint2_2e_t)(double*, int*, int*, int, int*, int, double*, void*);
typedef int (*cgto_t)(int, int*);

int  int1e_ipipovlp_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int1e_ipipovlp_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint1e_ipipovlp_cart(double*, int*, int*, int, int*, int, double*);
int  cint1e_ipipovlp_sph(double*, int*, int*, int, int*, int, double*);
int  int1e_ipovlpip_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int1e_ipovlpip_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint1e_ipovlpip_cart(double*, int*, int*, int, int*, int, double*);
int  cint1e_ipovlpip_sph(double*, int*, int*, int, int*, int, double*);
int  int1e_ipipkin_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int1e_ipipkin_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint1e_ipipkin_cart(double*, int*, int*, int, int*, int, double*);
int  cint1e_ipipkin_sph(double*, int*, int*, int, int*, int, double*);
int  int1e_ipkinip_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int1e_ipkinip_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint1e_ipkinip_cart(double*, int*, int*, int, int*, int, double*);
int  cint1e_ipkinip_sph(double*, int*, int*, int, int*, int, double*);
int  int1e_ipipnuc_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int1e_ipipnuc_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint1e_ipipnuc_cart(double*, int*, int*, int, int*, int, double*);
int  cint1e_ipipnuc_sph(double*, int*, int*, int, int*, int, double*);
int  int1e_ipnucip_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int1e_ipnucip_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint1e_ipnucip_cart(double*, int*, int*, int, int*, int, double*);
int  cint1e_ipnucip_sph(double*, int*, int*, int, int*, int, double*);
int  int2e_ipip1_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int2e_ipip1_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint2e_ipip1_cart(double*, int*, int*, int, int*, int, double*, void*);
int  cint2e_ipip1_sph(double*, int*, int*, int, int*, int, double*, void*);
int  int2e_ipvip1_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int2e_ipvip1_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint2e_ipvip1_cart(double*, int*, int*, int, int*, int, double*, void*);
int  cint2e_ipvip1_sph(double*, int*, int*, int, int*, int, double*, void*);
int  int2e_ip1ip2_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int2e_ip1ip2_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint2e_ip1ip2_cart(double*, int*, int*, int, int*, int, double*, void*);
int  cint2e_ip1ip2_sph(double*, int*, int*, int, int*, int, double*, void*);
int  int2e_ipip1ipip2_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int2e_ipip1ipip2_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint2e_ipip1ipip2_cart(double*, int*, int*, int, int*, int, double*, void*);
int  cint2e_ipip1ipip2_sph(double*, int*, int*, int, int*, int, double*, void*);
int  int2e_ipvip1ipvip2_cart(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  int2e_ipvip1ipvip2_sph(double*, int*, int*, int*, int, int*, int, double*, void*, double*);
int  cint2e_ipvip1ipvip2_cart(double*, int*, int*, int, int*, int, double*, void*);
int  cint2e_ipvip1ipvip2_sph(double*, int*, int*, int, int*, int, double*, void*);

/* Recorded from upstream libcint, same source, %.17g. */
#define REF_OVLP  2.0030499778307118
#define REF_KIN   5.5522155911991531
#define REF_NUC  -16.016935526669709
#define REF_E2    3.1717214410137018
#define REF_GMAX  0.052214489270357038
#define TOL       1e-14

#define NATM 3
#define NBAS 5
static int atm[NATM * 6], bas[NBAS * 8];
static double env[128];

static int fail;
static void expect(const char *what, double got, double want, double tol)
{
    if (!(fabs(got - want) <= tol)) {
        printf("  FAIL %-10s got %.17g want %.17g  (delta %.3g)\n",
               what, got, want, fabs(got - want));
        fail = 1;
    } else {
        printf("  ok   %-10s %.17g\n", what, got);
    }
}

/* Every entry above, walked the way the loops in main walk int2e: over every
   shell tuple the four-shell basis makes, comparing the CINT2 spelling with
   the CINT3 one element by element -- bitwise, since they are two names for
   one computation -- and reducing the block to two numbers that can be held
   against upstream.

   Two numbers rather than one because a derivative integral summed over the
   full index range can vanish by symmetry, exactly as int2e_ip1's does above,
   and a sum of zero is also what a routine that computed nothing would give.
   The largest element is what separates those two.  The table is 22 rows and
   one loop rather than 22 copies of the loop, and the reference pair on each
   row is what this same program printed built against upstream libcint. */
struct hess_case {
    const char *name;
    int         nshl, ncomp;
    cint3_t     f3;
    cint2_1e_t  f2_1e;   /* the CINT2 spelling of a one-electron integral */
    cint2_2e_t  f2_2e;   /* ... and of a two-electron one, which carries opt */
    cgto_t      dim;
    double      ref_sum, ref_max;
};

static const struct hess_case hess_cases[] = {
    {"int1e_ipipovlp_cart", 2, 9,
     int1e_ipipovlp_cart, cint1e_ipipovlp_cart, NULL,
     CINTcgto_cart, -21.764488341097277, 1.5046410473888836},
    {"int1e_ipipovlp_sph", 2, 9,
     int1e_ipipovlp_sph, cint1e_ipipovlp_sph, NULL,
     CINTcgto_spheric, -16.236883171551746, 0.94077916534999373},
    {"int1e_ipovlpip_cart", 2, 9,
     int1e_ipovlpip_cart, cint1e_ipovlpip_cart, NULL,
     CINTcgto_cart, 21.764488341097277, 1.5046410473888838},
    {"int1e_ipovlpip_sph", 2, 9,
     int1e_ipovlpip_sph, cint1e_ipovlpip_sph, NULL,
     CINTcgto_spheric, 16.236883171551746, 0.94077916534999395},
    {"int1e_ipipkin_cart", 2, 9,
     int1e_ipipkin_cart, cint1e_ipipkin_cart, NULL,
     CINTcgto_cart, -88.231406004752685, 6.1045436779777562},
    {"int1e_ipipkin_sph", 2, 9,
     int1e_ipipkin_sph, cint1e_ipipkin_sph, NULL,
     CINTcgto_spheric, -64.204421127991225, 3.3868049952599772},
    {"int1e_ipkinip_cart", 2, 9,
     int1e_ipkinip_cart, cint1e_ipkinip_cart, NULL,
     CINTcgto_cart, 88.231406004752657, 6.1045436779777589},
    {"int1e_ipkinip_sph", 2, 9,
     int1e_ipkinip_sph, cint1e_ipkinip_sph, NULL,
     CINTcgto_spheric, 64.204421127991225, 3.386804995259979},
    {"int1e_ipipnuc_cart", 2, 9,
     int1e_ipipnuc_cart, cint1e_ipipnuc_cart, NULL,
     CINTcgto_cart, 178.98647303354113, 11.665786955959351},
    {"int1e_ipipnuc_sph", 2, 9,
     int1e_ipipnuc_sph, cint1e_ipipnuc_sph, NULL,
     CINTcgto_spheric, 148.45039986927108, 7.8106279630790842},
    {"int1e_ipnucip_cart", 2, 9,
     int1e_ipnucip_cart, cint1e_ipnucip_cart, NULL,
     CINTcgto_cart, -198.27922836389132, 13.625461031231753},
    {"int1e_ipnucip_sph", 2, 9,
     int1e_ipnucip_sph, cint1e_ipnucip_sph, NULL,
     CINTcgto_spheric, -154.76782405917055, 8.1524815099052859},
    {"int2e_ipip1_cart", 4, 9,
     int2e_ipip1_cart, NULL, cint2e_ipip1_cart,
     CINTcgto_cart, -108.58376474596938, 1.0041936476453561},
    {"int2e_ipip1_sph", 4, 9,
     int2e_ipip1_sph, NULL, cint2e_ipip1_sph,
     CINTcgto_spheric, -26.459239108996645, 0.21659434322308249},
    {"int2e_ipvip1_cart", 4, 9,
     int2e_ipvip1_cart, NULL, cint2e_ipvip1_cart,
     CINTcgto_cart, 89.27619165259776, 0.85443300782272802},
    {"int2e_ipvip1_sph", 4, 9,
     int2e_ipvip1_sph, NULL, cint2e_ipvip1_sph,
     CINTcgto_spheric, 23.73663165553512, 0.20094146904507848},
    {"int2e_ip1ip2_cart", 4, 9,
     int2e_ip1ip2_cart, NULL, cint2e_ip1ip2_cart,
     CINTcgto_cart, 9.6537865466861721, 0.09141420238957311},
    {"int2e_ip1ip2_sph", 4, 9,
     int2e_ip1ip2_sph, NULL, cint2e_ip1ip2_sph,
     CINTcgto_spheric, 1.3613037267306964, 0.047829332609169821},
    {"int2e_ipip1ipip2_cart", 4, 81,
     int2e_ipip1ipip2_cart, NULL, cint2e_ipip1ipip2_cart,
     CINTcgto_cart, 534.37642842877551, 2.1303082892895389},
    {"int2e_ipip1ipip2_sph", 4, 81,
     int2e_ipip1ipip2_sph, NULL, cint2e_ipip1ipip2_sph,
     CINTcgto_spheric, 231.52855952630773, 0.68731789682667188},
    {"int2e_ipvip1ipvip2_cart", 4, 81,
     int2e_ipvip1ipvip2_cart, NULL, cint2e_ipvip1ipvip2_cart,
     CINTcgto_cart, 356.1177654216159, 1.669973712251964},
    {"int2e_ipvip1ipvip2_sph", 4, 81,
     int2e_ipvip1ipvip2_sph, NULL, cint2e_ipvip1ipvip2_sph,
     CINTcgto_spheric, 206.19955190471785, 0.63388989900075354},
};

/* 81 components over a d|d|d|d quartet in the Cartesian convention is
   81 * 6^4 doubles, and that is the largest block anything below asks for.
   The gradient blocks above fit in a fortieth of it. */
#define HBUF (81 * 6 * 6 * 6 * 6)
static double hbuf[HBUF], hbuf2[HBUF];

static void walk_hess(const struct hess_case *c)
{
    int idx, i, m, n, t, total = 1, shls[4], ok3, ok2;
    double sum = 0, mx = 0;
    char what[40];

    for (i = 0; i < c->nshl; i++) total *= NBAS;

    for (idx = 0; idx < total; idx++) {
        t = idx;
        n = c->ncomp;
        for (i = 0; i < c->nshl; i++) {
            shls[i] = t % NBAS;
            t /= NBAS;
            n *= c->dim(shls[i], bas);
        }
        if (n > HBUF) {
            printf("  FAIL %s: a block of %d does not fit the buffer\n",
                   c->name, n);
            fail = 1;
            return;
        }

        ok3 = c->f3(hbuf, NULL, shls, atm, NATM, bas, NBAS, env, NULL, NULL);
        ok2 = c->f2_1e ? c->f2_1e(hbuf2, shls, atm, NATM, bas, NBAS, env)
                       : c->f2_2e(hbuf2, shls, atm, NATM, bas, NBAS, env, NULL);
        if (ok3 != ok2) {
            printf("  FAIL %s: CINT2/CINT3 disagree on emptiness\n", c->name);
            fail = 1;
        }
        if (!ok3)
            continue;

        for (m = 0; m < n; m++) {
            if (hbuf[m] != hbuf2[m]) {
                printf("  FAIL %s: CINT2 != CINT3 at %d: %.17g vs %.17g\n",
                       c->name, m, hbuf[m], hbuf2[m]);
                fail = 1;
                break;
            }
            sum += hbuf[m];
            if (fabs(hbuf[m]) > mx) mx = fabs(hbuf[m]);
        }
    }

    /* TOL scaled by the magnitude, not applied flat.  A sum over Cartesian d
       quartets runs to several hundred, whose last bit is already 1e-13, and
       an absolute 1e-14 on that would be a check on whether two compilers add
       in the same order rather than on whether libfint reproduces libcint.
       Relative 1e-14 leaves the same room the rows above get from an absolute
       1e-14 on numbers of order one, and it is nowhere near loose enough to
       hide a real difference: a last-bit disagreement spread over a block of
       this size accumulates to sqrt(N) ulp, which is orders larger. */
    snprintf(what, sizeof what, "%s", c->name);
    expect(what, sum, c->ref_sum, TOL * fmax(1.0, fabs(c->ref_sum)));
    snprintf(what, sizeof what, "%s max", c->name);
    expect(what, mx, c->ref_max, TOL * fmax(1.0, fabs(c->ref_max)));
}

int main(void)
{
    int off = 20, i, j, k, l, m, e0, c0, e1, c1, ep, cp, ed, cd, shls[4], n;
    double buf[4096], buf2[4096];
    double s_ovlp = 0, s_kin = 0, s_nuc = 0, s_e2 = 0, s_e2opt = 0, gmax = 0;
    void *opt = NULL;

    memset(atm, 0, sizeof atm);
    memset(bas, 0, sizeof bas);
    memset(env, 0, sizeof env);

    /* Water-shaped: O at the origin, two H out along x and y, in Bohr. */
    atm[0 * 6 + 0] = 8; atm[0 * 6 + 1] = off;
    env[off++] = 0.0;  env[off++] = 0.0;  env[off++] = 0.0;
    atm[1 * 6 + 0] = 1; atm[1 * 6 + 1] = off;
    env[off++] = 1.43; env[off++] = 0.0;  env[off++] = 0.0;
    atm[2 * 6 + 0] = 1; atm[2 * 6 + 1] = off;
    env[off++] = 0.0;  env[off++] = 1.43; env[off++] = 0.0;

    e0 = off; env[off++] = 130.7; env[off++] = 23.81; env[off++] = 6.443;
    c0 = off; env[off++] = 0.154; env[off++] = 0.535; env[off++] = 0.445;
    e1 = off; env[off++] = 3.425; c1 = off; env[off++] = 1.0;
    ep = off; env[off++] = 1.0;   cp = off; env[off++] = 1.0;
    ed = off; env[off++] = 0.8;   cd = off; env[off++] = 1.0;

    /* A contracted s on O, a p on O, one s on each H, and a d back on O.  The
       first four give s, p and a contraction depth above one, which exercises
       every branch here; the d is what tells a Cartesian entry point from a
       spherical one at all.  Below l = 2 the two conventions agree value for
       value, so an export that quietly ran its twin would pass a basis of s
       and p and nothing would say so. */
    bas[0 * 8 + 0] = 0; bas[0 * 8 + 1] = 0; bas[0 * 8 + 2] = 3; bas[0 * 8 + 3] = 1;
    bas[0 * 8 + 5] = e0; bas[0 * 8 + 6] = c0;
    bas[1 * 8 + 0] = 0; bas[1 * 8 + 1] = 1; bas[1 * 8 + 2] = 1; bas[1 * 8 + 3] = 1;
    bas[1 * 8 + 5] = ep; bas[1 * 8 + 6] = cp;
    bas[2 * 8 + 0] = 1; bas[2 * 8 + 1] = 0; bas[2 * 8 + 2] = 1; bas[2 * 8 + 3] = 1;
    bas[2 * 8 + 5] = e1; bas[2 * 8 + 6] = c1;
    bas[3 * 8 + 0] = 2; bas[3 * 8 + 1] = 0; bas[3 * 8 + 2] = 1; bas[3 * 8 + 3] = 1;
    bas[3 * 8 + 5] = e1; bas[3 * 8 + 6] = c1;
    bas[4 * 8 + 0] = 0; bas[4 * 8 + 1] = 2; bas[4 * 8 + 2] = 1; bas[4 * 8 + 3] = 1;
    bas[4 * 8 + 5] = ed; bas[4 * 8 + 6] = cd;

    printf("c_abi_check: the C entry points, driven from C\n");

    if (CINTtot_cgto_spheric(bas, NBAS) != 11) {
        printf("  FAIL nao: got %d want 11\n", CINTtot_cgto_spheric(bas, NBAS));
        fail = 1;
    } else {
        printf("  ok   nao        11\n");
    }
    if (CINTcgto_spheric(1, bas) != 3) {
        printf("  FAIL cgto(p): got %d want 3\n", CINTcgto_spheric(1, bas));
        fail = 1;
    } else {
        printf("  ok   cgto(p)    3\n");
    }

    for (i = 0; i < NBAS; i++) for (j = 0; j < NBAS; j++) {
        shls[0] = i; shls[1] = j;
        n = CINTcgto_spheric(i, bas) * CINTcgto_spheric(j, bas);
        if (int1e_ovlp_sph(buf, NULL, shls, atm, NATM, bas, NBAS, env, NULL, NULL))
            for (m = 0; m < n; m++) s_ovlp += buf[m];
        if (int1e_kin_sph(buf, NULL, shls, atm, NATM, bas, NBAS, env, NULL, NULL))
            for (m = 0; m < n; m++) s_kin += buf[m];
        if (int1e_nuc_sph(buf, NULL, shls, atm, NATM, bas, NBAS, env, NULL, NULL))
            for (m = 0; m < n; m++) s_nuc += buf[m];
    }

    cint2e_sph_optimizer(&opt, atm, NATM, bas, NBAS, env);
    for (i = 0; i < NBAS; i++) for (j = 0; j < NBAS; j++)
    for (k = 0; k < NBAS; k++) for (l = 0; l < NBAS; l++) {
        int a, b, c;
        shls[0] = i; shls[1] = j; shls[2] = k; shls[3] = l;
        n = CINTcgto_spheric(i, bas) * CINTcgto_spheric(j, bas)
          * CINTcgto_spheric(k, bas) * CINTcgto_spheric(l, bas);

        a = int2e_sph(buf, NULL, shls, atm, NATM, bas, NBAS, env, NULL, NULL);
        b = cint2e_sph(buf2, shls, atm, NATM, bas, NBAS, env, NULL);
        if (a != b) { printf("  FAIL CINT2/CINT3 disagree on emptiness\n"); fail = 1; }
        if (a) for (m = 0; m < n; m++) {
            if (buf[m] != buf2[m]) {
                printf("  FAIL CINT2 != CINT3 at %d: %.17g vs %.17g\n",
                       m, buf[m], buf2[m]);
                fail = 1; break;
            }
            s_e2 += buf[m];
        }

        c = int2e_sph(buf2, NULL, shls, atm, NATM, bas, NBAS, env, opt, NULL);
        if (a != c) { printf("  FAIL opt/no-opt disagree on emptiness\n"); fail = 1; }
        if (c) for (m = 0; m < n; m++) {
            if (buf[m] != buf2[m]) {
                printf("  FAIL opt != no-opt at %d: %.17g vs %.17g\n",
                       m, buf[m], buf2[m]);
                fail = 1; break;
            }
            s_e2opt += buf2[m];
        }
    }
    CINTdel_optimizer(&opt);

    for (i = 0; i < NBAS; i++) for (j = 0; j < NBAS; j++)
    for (k = 0; k < NBAS; k++) for (l = 0; l < NBAS; l++) {
        shls[0] = i; shls[1] = j; shls[2] = k; shls[3] = l;
        n = 3 * CINTcgto_spheric(i, bas) * CINTcgto_spheric(j, bas)
              * CINTcgto_spheric(k, bas) * CINTcgto_spheric(l, bas);
        if (int2e_ip1_sph(buf, NULL, shls, atm, NATM, bas, NBAS, env, NULL, NULL))
            for (m = 0; m < n; m++)
                if (fabs(buf[m]) > gmax) gmax = fabs(buf[m]);
    }

    expect("ovlp",  s_ovlp,  REF_OVLP, TOL);
    expect("kin",   s_kin,   REF_KIN,  TOL);
    expect("nuc",   s_nuc,   REF_NUC,  TOL);
    expect("int2e", s_e2,    REF_E2,   TOL);
    expect("int2e+opt", s_e2opt, REF_E2, TOL);
    /* The signed sum of ip1 over the full index range is zero by translational
       invariance, so it would also be zero if nothing were computed at all.
       The largest element is what distinguishes those two. */
    expect("ip1max", gmax,   REF_GMAX, TOL);

    /* The second derivatives, which is what an analytic Hessian is built from
       and what nothing above this line touches. */
    for (i = 0; i < (int)(sizeof hess_cases / sizeof hess_cases[0]); i++)
        walk_hess(&hess_cases[i]);

    printf(fail ? "c_abi_check: FAILED\n" : "c_abi_check: all checks passed\n");
    return fail;
}
