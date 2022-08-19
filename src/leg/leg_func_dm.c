/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(leg_func_dm)(unsigned long nmax, const REAL *r, const REAL *ri,
                        REAL *dm)
{
    dm[0] = PREC(0.0);
    for (unsigned long n = 1; n <= nmax; n++)
        dm[n] = r[2 * n + 3] * ri[2 * n + 2];


    return;
}
