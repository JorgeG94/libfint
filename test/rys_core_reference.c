/*
 * Reference data for the Rys root finder and the Schmidt path (D3).
 *
 * Two sections.  The first builds coefficient tables whose rows are monic
 * polynomials with known interlacing real roots -- the structure R_dsmit
 * produces -- and runs _CINT_polynomial_roots over them; interlacing matters
 * because the fallback path in that function depends on it.  The second
 * sweeps all three C Schmidt variants so each Fortran precision can be
 * compared against its own counterpart.
 *
 * The Fortran quad path is compared against CINTqrys_schmidt, not
 * CINTlrys_schmidt: both are binary128, whereas the C's long double is 80-bit
 * x87 and disagrees for that reason alone.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MXR 32
int _CINT_polynomial_roots(double *roots, double *cs, int nroots);
int CINTrys_schmidt(int,double,double,double*,double*);
int CINTlrys_schmidt(int,double,double,double*,double*);
int CINTqrys_schmidt(int,double,double,double*,double*);


/* coefficients of prod_{j<k}(x - r_j), ascending order, into c[0..k] */
static void monic_from_roots(const double *r, int k, double *c)
{
        memset(c, 0, sizeof(double)*(k+1));
        c[0] = 1.0;
        int deg = 0;
        for (int j = 0; j < k; j++) {
                for (int i = deg+1; i >= 1; i--) c[i] = c[i-1] - r[j]*c[i];
                c[0] = -r[j]*c[0];
                deg++;
        }
}

int main(int argc, char **argv)
{
        FILE *o = fopen(argv[1], "wb");
        double cs[MXR*MXR], roots[MXR], r[MXR], c[MXR];
        int32_t n32;

        for (int nroots = 1; nroots <= 13; nroots++) {
                int nroots1 = nroots + 1;
                /* Chebyshev-like nodes in (0,1), nested so that the order-k
                 * set interlaces the order-(k+1) set. */
                for (int scenario = 0; scenario < 3; scenario++) {
                        memset(cs, 0, sizeof(cs));
                        for (int order = 0; order <= nroots; order++) {
                                for (int j = 0; j < order; j++) {
                                        double t = (j + 0.5) / order;
                                        if (scenario == 0) r[j] = t;
                                        else if (scenario == 1) r[j] = t*t;          /* clustered low */
                                        else r[j] = 1.0 - (1.0-t)*(1.0-t);           /* clustered high */
                                }
                                monic_from_roots(r, order, c);
                                for (int i = 0; i <= order; i++) cs[order*nroots1+i] = c[i];
                        }
                        memset(roots, 0, sizeof(roots));
                        int err = _CINT_polynomial_roots(roots, cs, nroots);

                        n32 = nroots;      fwrite(&n32, 4, 1, o);
                        n32 = scenario;    fwrite(&n32, 4, 1, o);
                        n32 = err;         fwrite(&n32, 4, 1, o);
                        fwrite(cs, sizeof(double), nroots1*nroots1, o);
                        fwrite(roots, sizeof(double), nroots, o);
                }
        }
        fclose(o);

        /* ---- section 2: the three Schmidt variants ---- */
        {
        FILE *o2 = fopen(argv[2], "wb");
        double u[64], w[64];
        double xs[] = {1e-9,1e-4,0.01,0.5,1.,2.5,5.,9.,10.,11.,18.,22.,30.,50.,60.,100.,200.};
        double ls[] = {0.,0.1,0.25,0.4,0.5,0.6,0.75,0.9,0.95,0.99};
        for (int n = 1; n <= 13; n++)
        for (unsigned i = 0; i < sizeof(xs)/sizeof(*xs); i++)
        for (unsigned j = 0; j < sizeof(ls)/sizeof(*ls); j++) {
                int32_t a[3];
                for (int variant = 0; variant < 3; variant++) {
                        memset(u, 0, sizeof(u)); memset(w, 0, sizeof(w));
                        int e = variant == 0 ? CINTrys_schmidt(n, xs[i], ls[j], u, w)
                              : variant == 1 ? CINTlrys_schmidt(n, xs[i], ls[j], u, w)
                                             : CINTqrys_schmidt(n, xs[i], ls[j], u, w);
                        a[0] = n; a[1] = variant; a[2] = e;
                        fwrite(a, 4, 3, o2); fwrite(&xs[i], 8, 1, o2); fwrite(&ls[j], 8, 1, o2);
                        fwrite(u, 8, n, o2); fwrite(w, 8, n, o2);
                }
        }
        fclose(o2);
        }
        return 0;
}
