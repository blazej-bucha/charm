/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_dda)(const CHARM(shc) *shcs1,
                    const CHARM(shc) *shcs2,
                    unsigned long nmax,
                    REAL *dda,
                    CHARM(err) *err)
{
    CHARM(shc_ddv)(shcs1, shcs2, nmax, dda, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }



    for (unsigned long n = 0; n <= nmax; n++)
        dda[n] = SQRT(dda[n]);


    return;
}
