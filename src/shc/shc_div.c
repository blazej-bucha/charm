/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../misc/misc_is_nearly_equal.h"
#include "shc_arithmetics_checks.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_div)(CHARM(shc) *rop,
                    CHARM(shc) *op1,
                    CHARM(shc) *op2,
                    unsigned long nmin,
                    unsigned long nmax,
                    CHARM(err) *err)
{
    /* Sanity checks */
    /* --------------------------------------------------------------------- */
    CHARM(shc_arithmetics_checks)(rop, op1, op2, nmin, nmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    if (nmax > op2->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFUNCARG,
                       "\"nmax\" cannot be larger than \"op2->nmax\", as "
                       "this requires zero padding \"op2\" beyond "
                       "\"nmax\" which subsequently leads to division by "
                       "zero.");
        return;
    }


    /* Check that none of the coefficients in "op2" are zero within our degree
     * range from "nmin" to "nmax". */
    unsigned long nmm;
    _Bool czero, szero;
    for (unsigned long m = 0; m <= nmax; m++)
    {
        for (unsigned long n = CHARM_MAX(nmin, m); n <= nmax; n++)
        {
            nmm = n - m;


            czero = CHARM(misc_is_nearly_equal)(op2->c[m][nmm], PREC(0.0),
                                                CHARM(glob_threshold));
            szero = CHARM(misc_is_nearly_equal)(op2->s[m][nmm], PREC(0.0),
                                                CHARM(glob_threshold));
            if (czero || (szero && (m > 0)))
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG,
                               "The coefficients of \"op2\" must not be zero "
                               "in order to avoid division by zero.  The only "
                               "exceptions are the non-existing \"s\" "
                               "coefficients of order \"0\".");
                return;
            }
        }
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    for (unsigned long m = 0; m <= nmax; m++)
    {
        for (unsigned long n = CHARM_MAX(nmin, m); n <= nmax; n++)
        {
            nmm = n - m;


            if ((n <= op1->nmax) && (n <= op2->nmax))
            {
                rop->c[m][nmm] = op1->c[m][nmm] / op2->c[m][nmm];
                if (m > 0)
                    rop->s[m][nmm] = op1->s[m][nmm] / op2->s[m][nmm];
            }
            else
            {
                rop->c[m][nmm] = PREC(0.0);
                rop->s[m][nmm] = PREC(0.0);
            }
        }
    }
    /* --------------------------------------------------------------------- */


    return;
}
