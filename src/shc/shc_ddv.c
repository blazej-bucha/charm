/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_check_distribution.h"
#include "../misc/misc_is_nearly_equal.h"
#include "shc_check_distribution.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_ddv)(const CHARM(shc) *shcs1,
                    const CHARM(shc) *shcs2,
                    unsigned long nmax,
                    REAL *ddv,
                    CHARM(err) *err)
{
    /* Some error checks */
    /* --------------------------------------------------------------------- */
    CHARM(err_check_distribution)(err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(shc_check_distribution)(shcs1, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(shc_check_distribution)(shcs2, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    if (nmax > shcs1->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmax\" cannot be larger than \"shcs1->nmax\".");
        return;
    }


    if (nmax > shcs2->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmax\" cannot be larger than \"shcs2->nmax\".");
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(shcs1->r, shcs2->r,
                                     CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "To get meaningful difference degree "
                       "variances/amplitudes, \"shcs1->r\" must be equal "
                       "to \"shcs2->r\".");
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(shcs1->mu, shcs2->mu,
                                     CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "To get meaningful difference degree "
                       "variances/amplitudes, \"shcs1->mu\" must be equal "
                       "to \"shcs2->mu\".");
        return;
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    REAL ddvn, tmp_cnm, tmp_snm;
    for (unsigned long n = 0; n <= nmax; n++)
    {
        ddvn = PREC(0.0);
        for (unsigned long m = 0; m <= n; m++)
        {
            tmp_cnm = shcs1->c[m][n - m] - shcs2->c[m][n - m];
            tmp_snm = shcs1->s[m][n - m] - shcs2->s[m][n - m];
            ddvn += (tmp_cnm * tmp_cnm) + (tmp_snm * tmp_snm);
        }


        ddv[n] = ddvn;
    }
    /* --------------------------------------------------------------------- */


    return;
}
