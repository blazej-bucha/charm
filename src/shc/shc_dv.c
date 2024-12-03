/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_check_distribution.h"
#include "shc_check_distribution.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_dv)(const CHARM(shc) *shcs,
                   unsigned long nmax,
                   REAL *dv,
                   CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    CHARM(err_check_distribution)(err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(shc_check_distribution)(shcs, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Check the maximum harmonic degree */
    /* --------------------------------------------------------------------- */
    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmax\" cannot be larger than \"shcs->nmax\".");


        return;
    }
    /* --------------------------------------------------------------------- */


    REAL dvn;
    for (unsigned long n = 0; n <= nmax; n++)
    {
        dvn = PREC(0.0);
        for (unsigned long m = 0; m <= n; m++)
            dvn += shcs->c[m][n - m] * shcs->c[m][n - m] +
                   shcs->s[m][n - m] * shcs->s[m][n - m];


        dv[n] = dvn;
    }


    return;
}
