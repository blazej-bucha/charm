/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






/* An internal function to compute "POW[n * SIMD_SIZE + v] = (rref / r)^(n
 * + 1)" for "n = 0, 1, ..., nmax". */
void CHARM(shs_rpows)(size_t v, REAL rref, REAL r, REAL *pow,
                      unsigned long nmax)
{
    REAL ratio = rref / r;


    pow[v] = ratio;
    if (nmax == 0)
        return;


    REAL pow_tmp = ratio * ratio;
    for (unsigned long n = 1; n <= nmax; n++)
    {
        pow[n * SIMD_SIZE + v]  = pow_tmp;
        pow_tmp                *= ratio;
    }


    return;
}

