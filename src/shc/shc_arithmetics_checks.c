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
#include "shc_arithmetics_checks.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_arithmetics_checks)(CHARM(shc) *rop,
                                   CHARM(shc) *op1,
                                   CHARM(shc) *op2,
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


    CHARM(shc_check_distribution)(rop, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(shc_check_distribution)(op1, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(shc_check_distribution)(op2, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(rop->mu, op1->mu, CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"rop->mu\" must be equal to \"op1->mu\".");
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(rop->mu, op2->mu, CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"rop->mu\" must be equal to \"op2->mu\".");
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(rop->r, op1->r, CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"rop->r\" must be equal to \"op1->r\".");
        return;
    }


    if (!CHARM(misc_is_nearly_equal)(rop->r, op2->r, CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"rop->r\" must be equal to \"op2->r\".");
        return;
    }


    if (nmin > nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmin\" cannot be larger than \"nmax\".");
        return;
    }


    if (nmin > rop->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmin\" cannot be larger than \"rop->nmax\".");
        return;
    }


    if (nmax > rop->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmax\" cannot be larger than \"rop->nmax\".");
        return;
    }


    return;
}
