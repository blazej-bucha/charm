/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Computes sectorial Legendre functions for internal purposes.  The function
 * is based on the "alfsp" Fortran subroutine (Table 2) due to Fukushima
 * (2012). */
void CHARM(leg_func_prepare)(REAL u, REAL *ps, int *ips, const REAL *dm,
                             unsigned long nmax)
{
    /* No sectorial Legendre functions for maximum harmonic degree "0" */
    if (nmax == 0)
        return;


    REAL y;
    REAL x = ROOT3 * u;
    int ix = 0;


     ps[0] = x;
    ips[0] = ix;


    for (unsigned long n = 1; n < nmax; n++)
    {
        x = (dm[n] * u) * x;
        y = FABS(x);


        if (y >= BIGS)
        {
             x *= BIGI;
            ix += 1;
        }
        else if (y < BIGSI)
        {
             x *= BIG;
            ix -= 1;
        }


         ps[n] = x;
        ips[n] = ix;
    }


    return;
}
