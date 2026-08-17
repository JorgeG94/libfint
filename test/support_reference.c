/*
 * Reference data for the support layer (D4): the Cartesian-to-spherical
 * transforms, CINTgto_norm, and the small array helpers of fblas.c.
 *
 * The per-l transforms are static in the C, so they are reached through the
 * exported dispatch tables c2s_bra_sph[] / c2s_ket_sph[] rather than by name.
 * lds is deliberately set to nbra+2 so a routine that ignored the stride
 * would be caught.
 *
 * Every one of these should come out bit-identical: the transforms are the
 * same sparse sums in the same order, and a transpose is a permutation.
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
extern double *(*c2s_bra_sph[])(double *gsph, int nket, double *gcart, int l);
extern double *(*c2s_ket_sph[])(double *gsph, double *gcart, int lds, int nbra, int l);
double CINTgto_norm(int n, double a);
int CINTcgto_cart(int bas_id, const int *bas);
void CINTdmat_transpose(double*,double*,int,int);
void CINTdplus_transpose(double*,double*,int,int);
void CINTdaxpy2v(int,double,double*,double*,double*);


int main(int argc, char **argv)
{
        FILE *o = fopen(argv[1], "wb");
        double gcart[4096], gsph[4096], *r;
        int32_t hdr[3];

        for (int l = 0; l <= 6; l++) {
                int nf = (l+1)*(l+2)/2;
                int nd = 2*l+1;
                for (int nket = 1; nket <= 4; nket++) {
                        /* deterministic, spread over several magnitudes so a
                         * dropped term cannot hide */
                        for (int i = 0; i < nf*nket; i++) {
                                gcart[i] = (i % 7 + 1) * 0.37 * ((i % 3) ? 1 : -1)
                                         * (1.0 + 0.01*i);
                        }
                        memset(gsph, 0, sizeof(gsph));
                        r = c2s_bra_sph[l](gsph, nket, gcart, l);
                        hdr[0] = 0; hdr[1] = l; hdr[2] = nket;
                        fwrite(hdr, 4, 3, o);
                        fwrite(gcart, sizeof(double), nf*nket, o);
                        fwrite(r, sizeof(double), nd*nket, o);
                }
                for (int nbra = 1; nbra <= 4; nbra++) {
                        int lds = nbra + 2;   /* deliberately != nbra */
                        for (int i = 0; i < nf*nbra; i++) {
                                gcart[i] = (i % 5 + 1) * 0.53 * ((i % 4) ? 1 : -1)
                                         * (1.0 + 0.02*i);
                        }
                        memset(gsph, 0, sizeof(gsph));
                        r = c2s_ket_sph[l](gsph, gcart, lds, nbra, l);
                        hdr[0] = 1; hdr[1] = l; hdr[2] = nbra;
                        fwrite(hdr, 4, 3, o);
                        fwrite(gcart, sizeof(double), nf*nbra, o);
                        /* the C may return gcart; read nd*lds from wherever it points */
                        fwrite(r, sizeof(double), nd*lds, o);
                }
        }
        /* gto_norm over a spread of exponents and l */
        for (int l = 0; l <= 8; l++) {
                for (int k = 0; k < 12; k++) {
                        double a = pow(10.0, k - 5);
                        double v = CINTgto_norm(l, a);
                        hdr[0] = 2; hdr[1] = l; hdr[2] = k;
                        fwrite(hdr, 4, 3, o);
                        fwrite(&a, sizeof(double), 1, o);
                        fwrite(&v, sizeof(double), 1, o);
                }
        }
        fclose(o);

        /* ---- the fblas helpers ---- */
        {
        FILE *o2 = fopen(argv[2], "wb");
        double aa[4096], tt[4096], yy[4096], vv[4096];
        int32_t h[3];
        for (int m = 1; m <= 17; m++) for (int n = 1; n <= 17; n++) {
                for (int i = 0; i < m*n; i++) {
                        aa[i] = (i%9+1)*0.31*((i%2)?1:-1);
                        tt[i] = (i%5+1)*0.11;
                }
                CINTdmat_transpose(tt, aa, m, n);
                h[0]=0; h[1]=m; h[2]=n;
                fwrite(h,4,3,o2); fwrite(aa,8,m*n,o2); fwrite(tt,8,m*n,o2);
                for (int i = 0; i < m*n; i++) tt[i] = (i%5+1)*0.11;
                CINTdplus_transpose(tt, aa, m, n);
                h[0]=1;
                fwrite(h,4,3,o2); fwrite(aa,8,m*n,o2); fwrite(tt,8,m*n,o2);
        }
        for (int n = 1; n <= 64; n++) {
                for (int i = 0; i < n; i++) {
                        aa[i] = (i%7+1)*0.29;
                        yy[i] = (i%3+1)*0.71;
                }
                CINTdaxpy2v(n, 1.37, aa, yy, vv);
                h[0]=2; h[1]=n; h[2]=0;
                fwrite(h,4,3,o2); fwrite(aa,8,n,o2); fwrite(yy,8,n,o2); fwrite(vv,8,n,o2);
        }
        fclose(o2);
        }
        return 0;
}
