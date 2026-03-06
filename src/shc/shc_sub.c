/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "shc_arithmetics_checks.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_sub)(CHARM(shc) *rop,
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
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    unsigned long nmm;
    _Bool ngtop1, ngtop2;
    for (unsigned long m = 0; m <= nmax; m++)
    {
        for (unsigned long n = CHARM_MAX(nmin, m); n <= nmax; n++)
        {
            nmm = n - m;
            ngtop1 = n > op1->nmax;
            ngtop2 = n > op2->nmax;


            if (ngtop1 && ngtop2)
            {
                rop->c[m][nmm] = PREC(0.0);
                rop->s[m][nmm] = PREC(0.0);
            }
            else if (ngtop1)
            {
                rop->c[m][nmm] = -op2->c[m][nmm];
                rop->s[m][nmm] = -op2->s[m][nmm];
            }
            else if (ngtop2)
            {
                rop->c[m][nmm] = op1->c[m][nmm];
                rop->s[m][nmm] = op1->s[m][nmm];
            }
            else
            {
                rop->c[m][nmm] = op1->c[m][nmm] - op2->c[m][nmm];
                rop->s[m][nmm] = op1->s[m][nmm] - op2->s[m][nmm];
            }
        }
    }
    /* --------------------------------------------------------------------- */


    return;
}
