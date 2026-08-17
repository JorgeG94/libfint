/*
 * Reference generator for the Rys quadrature port (deliverable D3).
 *
 * Dumps roots and weights from the C implementation across a sweep chosen to
 * straddle every branch constant in CINTrys_roots and CINTsr_rys_roots, so
 * that a Fortran reimplementation is checked on each arm of the ladder rather
 * than on whatever happens to be common.
 *
 * The branch structure being covered:
 *
 *   CINTrys_roots(nroots, x)
 *     x <= 3e-7                  SMALLX_LIMIT polynomial shortcut
 *     x >= 35 + nroots*5         large-x asymptotic shortcut
 *     nroots 1..5                closed forms rys_root1..5
 *     nroots 6,7                 jacobi / schmidt   at x = 11
 *     nroots 8                   jacobi / lschmidt  at x = 11
 *     nroots 9                   ljacobi / llaguerre at x = 10
 *     nroots 10,11               ljacobi / llaguerre at x = 18
 *     nroots 12                  ljacobi / llaguerre at x = 22
 *     nroots >12                 qjacobi / qlaguerre at x = 50
 *
 *   CINTsr_rys_roots(nroots, x, lower) adds breakpoints on `lower` at
 *   0.15, 0.25, 0.4, 0.5, 0.6, 0.65, 0.75, 0.8, 0.9, 0.93, 0.97, 0.99.
 *
 * A note on how far the sweep reaches.  CINTsr_rys_roots calls exit() on a
 * failed solve unless libcint itself was built with -DKEEP_GOING -- and the
 * exit is inside the library, so defining it here would not help.  The sweep
 * therefore stops at x = 200 by default, which is below where the C starts
 * failing.  Pass a third argument to raise the cap when the library has been
 * built with KEEP_GOING, and the sweep will then also cover the region where
 * the C gives up.
 *
 * Output is a binary stream of records, little-endian:
 *     int32 kind      0 = long range, 1 = short range
 *     int32 nroots
 *     double x
 *     double lower    (0 for kind 0)
 *     int32 err
 *     double u[nroots]
 *     double w[nroots]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cint_config.h"
#include "rys_roots.h"

#define MAXN 32

static FILE *out;

static void emit(int kind, int nroots, double x, double lower, int err,
                 const double *u, const double *w)
{
        int32_t k = kind, n = nroots, e = err;
        fwrite(&k, sizeof(int32_t), 1, out);
        fwrite(&n, sizeof(int32_t), 1, out);
        fwrite(&x, sizeof(double), 1, out);
        fwrite(&lower, sizeof(double), 1, out);
        fwrite(&e, sizeof(int32_t), 1, out);
        fwrite(u, sizeof(double), nroots, out);
        fwrite(w, sizeof(double), nroots, out);
}

/* x values that bracket every shortcut and segment breakpoint for a given
 * nroots.  Straddling matters more than density: a value just below and just
 * above each constant exercises both arms. */
static double g_xmax = 200.;

static int build_x_list(int nroots, double *xs)
{
        static const double bp[] = {10., 11., 18., 22., 50., 60.};
        int n = 0;
        double largex = 35. + nroots * 5.;

        xs[n++] = 0.0;
        xs[n++] = 1e-9;             /* below SMALLX_LIMIT */
        xs[n++] = 3e-7;             /* exactly SMALLX_LIMIT */
        xs[n++] = 3.1e-7;           /* just above */
        xs[n++] = 1e-4;
        xs[n++] = 0.01;
        xs[n++] = 0.5;
        xs[n++] = 1.0;
        xs[n++] = 2.5;
        xs[n++] = 5.0;
        for (unsigned i = 0; i < sizeof(bp) / sizeof(bp[0]); i++) {
                xs[n++] = bp[i] - 1e-6;
                xs[n++] = bp[i];
                xs[n++] = bp[i] + 1e-6;
        }
        xs[n++] = largex - 1e-6;    /* just inside the solver */
        xs[n++] = largex;           /* exactly the asymptotic cut */
        xs[n++] = largex + 1e-6;    /* just into the asymptotic branch */
        xs[n++] = largex * 2;
        xs[n++] = 200.;
        xs[n++] = 1000.;

        /* drop anything past the cap, keeping the order stable */
        int m = 0;
        for (int i = 0; i < n; i++) {
                if (xs[i] <= g_xmax) {
                        xs[m++] = xs[i];
                }
        }
        return m;
}

int main(int argc, char **argv)
{
        double u[MAXN], w[MAXN], xs[64];
        static const double lowers[] = {
                0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.55, 0.6, 0.65,
                0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 0.995,
        };
        int max_nroots = 14;
        long nrec = 0;

        if (argc < 2) {
                fprintf(stderr, "usage: rys_reference <outfile> [max_nroots] [max_x]\n");
                return 2;
        }
        if (argc > 2) {
                max_nroots = atoi(argv[2]);
        }
        if (argc > 3) {
                g_xmax = atof(argv[3]);
        }
        out = fopen(argv[1], "wb");
        if (!out) {
                fprintf(stderr, "cannot open %s\n", argv[1]);
                return 2;
        }

        for (int nroots = 1; nroots <= max_nroots; nroots++) {
                int nx = build_x_list(nroots, xs);

                for (int i = 0; i < nx; i++) {
                        memset(u, 0, sizeof(u));
                        memset(w, 0, sizeof(w));
                        int err = CINTrys_roots(nroots, xs[i], u, w);
                        emit(0, nroots, xs[i], 0.0, err, u, w);
                        nrec++;
                }

                for (int i = 0; i < nx; i++) {
                        for (unsigned j = 0; j < sizeof(lowers) / sizeof(lowers[0]); j++) {
                                memset(u, 0, sizeof(u));
                                memset(w, 0, sizeof(w));
                                int err = CINTsr_rys_roots(nroots, xs[i], lowers[j], u, w);
                                emit(1, nroots, xs[i], lowers[j], err, u, w);
                                nrec++;
                        }
                }
        }

        fclose(out);
        fprintf(stderr, "wrote %ld records, nroots 1..%d, x <= %g\n",
                nrec, max_nroots, g_xmax);
        return 0;
}
