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
int  CINTtot_cgto_spheric(int*, int);

/* Recorded from upstream libcint, same source, %.17g. */
#define REF_OVLP  0.48086061770996963
#define REF_KIN   1.2916778339954504
#define REF_NUC  -4.8028584141735244
#define REF_E2    0.20182398154332476
#define REF_GMAX  0.003606113676330098
#define TOL       1e-14

#define NATM 3
#define NBAS 4
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

int main(void)
{
    int off = 20, i, j, k, l, m, e0, c0, e1, c1, ep, cp, shls[4], n;
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

    /* A contracted s on O, a p on O, one s on each H: s, p and a contraction
       depth above one, which is enough to exercise every branch here. */
    bas[0 * 8 + 0] = 0; bas[0 * 8 + 1] = 0; bas[0 * 8 + 2] = 3; bas[0 * 8 + 3] = 1;
    bas[0 * 8 + 5] = e0; bas[0 * 8 + 6] = c0;
    bas[1 * 8 + 0] = 0; bas[1 * 8 + 1] = 1; bas[1 * 8 + 2] = 1; bas[1 * 8 + 3] = 1;
    bas[1 * 8 + 5] = ep; bas[1 * 8 + 6] = cp;
    bas[2 * 8 + 0] = 1; bas[2 * 8 + 1] = 0; bas[2 * 8 + 2] = 1; bas[2 * 8 + 3] = 1;
    bas[2 * 8 + 5] = e1; bas[2 * 8 + 6] = c1;
    bas[3 * 8 + 0] = 2; bas[3 * 8 + 1] = 0; bas[3 * 8 + 2] = 1; bas[3 * 8 + 3] = 1;
    bas[3 * 8 + 5] = e1; bas[3 * 8 + 6] = c1;

    printf("c_abi_check: the C entry points, driven from C\n");

    if (CINTtot_cgto_spheric(bas, NBAS) != 6) {
        printf("  FAIL nao: got %d want 6\n", CINTtot_cgto_spheric(bas, NBAS));
        fail = 1;
    } else {
        printf("  ok   nao        6\n");
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

    printf(fail ? "c_abi_check: FAILED\n" : "c_abi_check: all checks passed\n");
    return fail;
}
