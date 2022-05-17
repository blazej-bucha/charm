/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* An internal function to compute "POW[n] = (rref / r)^(n + 1)" for "n = 0, 1,
 * ..., nmax". */
void CHARM(shs_rpows)(REAL rref, REAL r, REAL *pow, unsigned long nmax)
{
    REAL ratio   = rref  / r;


    pow[0] = ratio;
    REAL pow_tmp = ratio * ratio;
    for (unsigned long n = 1; n <= nmax; n++)
    {
        pow[n]   = pow_tmp;
        pow_tmp *= ratio;
    }


    return;
}
