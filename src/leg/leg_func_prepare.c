/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






/* Computes sectorial Legendre functions for internal purposes.  The function
 * is based on the "alfsp" Fortran subroutine (Table 2) due to Fukushima
 * (2012). */
void CHARM(leg_func_prepare)(const REAL *u, REAL *ps, int *ips, const REAL *dm,
                             unsigned long nmax)
{
    /* No sectorial Legendre functions for maximum harmonic degree "0" */
    if (nmax == 0)
        return;


    REAL y[SIMD_SIZE], x[SIMD_SIZE];
    int ix[SIMD_SIZE];
    for (size_t v = 0; v < SIMD_SIZE; v++)
    {
         x[v] = ROOT3 * u[v];
        ix[v] = 0;


         ps[v] =  x[v];
        ips[v] = ix[v];
    }


    for (unsigned long n = 1; n < nmax; n++)
    {
        for (size_t v = 0; v < SIMD_SIZE; v++)
        {
            x[v] = (dm[n] * u[v]) * x[v];
            y[v] = FABS(x[v]);


            if (y[v] >= BIGS)
            {
                 x[v] *= BIGI;
                ix[v] += 1;
            }
            else if (y[v] < BIGSI)
            {
                 x[v] *= BIG;
                ix[v] -= 1;
            }


             ps[n * SIMD_SIZE + v] =  x[v];
            ips[n * SIMD_SIZE + v] = ix[v];
        }
    }


    return;
}

