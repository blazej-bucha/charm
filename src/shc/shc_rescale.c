/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_rescale)(CHARM(shc) *shcs, REAL munew, REAL rnew,
           CHARM(err) *err)
{
    /* Pre-computations */
    /* --------------------------------------------------------------------- */
    REAL mu_ratio = shcs->mu / munew;
    REAL r_ratio  = shcs->r  / rnew;


    REAL *tmp = (REAL *)malloc((shcs->nmax + 1) * sizeof(REAL));
    if (tmp == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        return;
    }


    REAL tmp2 = mu_ratio * r_ratio;
    tmp[0] = mu_ratio;
    for (unsigned long n = 1; n <= shcs->nmax; n++)
    {
        tmp[n] = tmp2;
        tmp2  *= r_ratio;
    }
    /* --------------------------------------------------------------------- */






    /* Rescale the coefficients */
    /* --------------------------------------------------------------------- */
    for (unsigned long m = 0; m <= shcs->nmax; m++)
    {
        for (unsigned long n = m; n <= shcs->nmax; n++)
        {
            shcs->c[m][n - m] *= tmp[n];
            shcs->s[m][n - m] *= tmp[n];
        }
    }


    free(tmp);
    /* --------------------------------------------------------------------- */






    /* Finally, update the "mu" and "r" members of "shcs" to the new values */
    /* --------------------------------------------------------------------- */
    shcs->mu = munew;
    shcs->r  = rnew;
    /* --------------------------------------------------------------------- */






    return;
}
