/*
 * Dumps the C's own Boys function so the Fortran double path can be checked
 * against it value by value (deliverable D3).
 *
 * Only the double routines are dumped.  The C's long double path is not a
 * useful reference for the Fortran quadruple one: at t = 1e-12, lower = 0.75,
 * m = 1 the exact F_1 is 0.19270833..., the C double gives -11.12 and the C
 * long double 0.2142, against 0.1927083333332 from the Fortran.  The quad path
 * is therefore checked against the Boys recurrence instead, in fmt_check.F90.
 *
 * Layout, little-endian, for m = 0..20 and each t in the sweep:
 *   double gamma_inc_like[m+1]
 *   then for each lower: double fmt_erfc_like[m+1]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cint_config.h"
#include "rys_roots.h"
int main(int argc, char**argv){
    FILE*o=fopen(argv[1],"wb");
    double f[64];
    /* t values straddling the turnover points and both asymptotic regimes */
    double ts[]={0.,1e-12,1e-6,1e-3,0.01,0.1,0.5,0.86602540378,0.87,1.0,1.295010032056,
                 1.3,2.0,3.28,5.0,8.0,12.0,20.0,40.0,100.0,200.0,500.0};
    double ls[]={0.,0.1,0.25,0.5,0.75,0.9,0.95,0.99};
    for(int m=0;m<=20;m++)
    for(unsigned i=0;i<sizeof(ts)/sizeof(*ts);i++){
        memset(f,0,sizeof(f));
        gamma_inc_like(f,ts[i],m);
        fwrite(f,sizeof(double),m+1,o);
        for(unsigned j=0;j<sizeof(ls)/sizeof(*ls);j++){
            memset(f,0,sizeof(f));
            fmt_erfc_like(f,ts[i],ls[j],m);
            fwrite(f,sizeof(double),m+1,o);
        }
    }
    fclose(o); return 0;
}
