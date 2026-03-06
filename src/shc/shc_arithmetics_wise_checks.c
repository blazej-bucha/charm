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
#include "shc_arithmetics_wise_checks.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_arithmetics_wise_checks)(CHARM(shc) *shcs,
                                        unsigned long nmin,
                                        unsigned long nmax,
                                        CHARM(err) *err)
{
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


    if (nmin > nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmin\" cannot be larger than \"nmax\".");
        return;
    }


    if (nmin > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmin\" cannot be larger than \"shcs->nmax\".");
        return;
    }


    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmax\" cannot be larger than \"shcs->nmax\".");
        return;
    }


    return;
}
