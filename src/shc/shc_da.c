/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_da)(const CHARM(shc) *shcs, unsigned long nmax, REAL *da,
                   CHARM(err) *err)
{
    CHARM(shc_dv)(shcs, nmax, da, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }



    for (unsigned long n = 0; n <= nmax; n++)
        da[n] = SQRT(da[n]);


    return;
}
