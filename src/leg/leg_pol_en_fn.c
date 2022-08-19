/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(leg_pol_en_fn)(unsigned long nmax, REAL *en, REAL *fn)
{
    en[0] = fn[0] = PREC(0.0);


    if (nmax == 0)
        return;


    en[0] = en[1] = fn[0] = fn[1] = PREC(0.0);
    if (nmax == 1)
        return;


    REAL n_flt;
    for (unsigned long n = 2; n <= nmax; n++)
    {
        n_flt = (REAL)n;
        en[n] = (REAL)(2 * n - 1) / n_flt;
        fn[n] = (REAL)(n - 1) / n_flt;
    }


    return;
}
