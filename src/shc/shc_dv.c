/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_dv)(const CHARM(shc) *shcs, unsigned long nmax, REAL *dv,
                   CHARM(err) *err)
{
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
