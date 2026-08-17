/* Dump _CINTdiagonalize's answers for a spread of tridiagonal matrices, so a
 * Fortran translation can be diffed against them bit for bit. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define MXRYSROOTS 32
int _CINTdiagonalize(int n, double *diag, double *diag_off1, double *eig, double *vec);

/* deterministic, portable: no rand() so the two sides see identical input */
static unsigned long long st;
static double nextd(void) {
    st = st * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)((st >> 11) & ((1ULL << 53) - 1)) / (double)(1ULL << 53);
}

int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "wb");
    double d[MXRYSROOTS], e[MXRYSROOTS], w[MXRYSROOTS], z[MXRYSROOTS*MXRYSROOTS];
    int n, k, i, info, ncase = 0;
    st = 88172645463325252ULL;
    for (n = 1; n <= 16; n++) {
      for (k = 0; k < 40; k++) {
        /* Jacobi-like: positive off-diagonals, which is what Wheeler feeds it */
        for (i = 0; i < n; i++) d[i] = 0.5 * nextd() + (double)i * 0.25;
        for (i = 0; i < n; i++) e[i] = 0.1 + 0.9 * nextd();
        if (k % 4 == 1) for (i = 0; i < n; i++) e[i] *= 1e-6;   /* near-decoupled */
        if (k % 4 == 2) for (i = 0; i < n; i++) d[i] = 1.0;     /* clustered */
        if (k % 4 == 3) for (i = 0; i < n; i++) e[i] *= 1e3;    /* dominant */
        memset(w, 0, sizeof w); memset(z, 0, sizeof z);
        /* the INPUT goes to the file too, so the Fortran side never has to
         * reproduce this generator -- only read what it produced */
        fwrite(&n, sizeof n, 1, f);
        fwrite(d, sizeof(double), n, f);
        fwrite(e, sizeof(double), n, f);
        info = _CINTdiagonalize(n, d, e, w, z);
        fwrite(&info, sizeof info, 1, f);
        fwrite(w, sizeof(double), n, f);
        fwrite(z, sizeof(double), (size_t)n*n, f);
        ncase++;
      }
    }
    fclose(f);
    printf("  %d cases -> %s\n", ncase, argv[1]);
    return 0;
}
