/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "shc_arithmetics_wise_checks.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_mul_degree_wise)(CHARM(shc) *shcs,
                                const REAL *a,
                                unsigned long nmin,
                                unsigned long nmax,
                                CHARM(err) *err)
{
    /* Sanity checks */
    /* --------------------------------------------------------------------- */
    CHARM(shc_arithmetics_wise_checks)(shcs, nmin, nmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    unsigned long nmm, i;
    for (unsigned long m = 0; m <= nmax; m++)
    {
        for (unsigned long n = CHARM_MAX(nmin, m); n <= nmax; n++)
        {
            nmm = n - m;
            i = n - nmin;
            shcs->c[m][nmm] *= a[i];
            shcs->s[m][nmm] *= a[i];
        }
    }
    /* --------------------------------------------------------------------- */


    return;
}
