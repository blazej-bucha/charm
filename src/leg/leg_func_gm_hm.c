/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(leg_func_gm_hm)(unsigned long nmax, const REAL *r,
                           const REAL *ri, REAL *gm, REAL *hm)
{
    gm[0] = ADDP(0.0);
    if (nmax == 0)
        return;


    gm[1] = ADDP(0.0);
    if (nmax == 1)
        return;


    gm[2] = ADDP(0.0);
    if (nmax == 2)
        return;


    for (unsigned long n = 3; n <= nmax; n++)
        gm[n] = ADDP(1.0) / (REAL)(2 * n + 2) * r[n] * r[2 * n + 1]
              * r[2 * n - 1] * ri[n - 1];


    for (unsigned long n = 0; n <= nmax; n++)
        hm[n] = (REAL)(n - 2) / (REAL)(n + 1);


    return;
}
