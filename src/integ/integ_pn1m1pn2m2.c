/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "integ_ccs.h"
#include "integ_css.h"
#include "integ_scs.h"
#include "integ_sss.h"
#include "../leg/leg_pnmj_check_ordering.h"
#include "../misc/misc_nan.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_check_distribution.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(integ_pn1m1pn2m2)(REAL cltmin,
                             REAL cltmax,
                             unsigned long n1,
                             unsigned long m1,
                             unsigned long n2,
                             unsigned long m2,
                             const CHARM(pnmj) *pnmj,
                             CHARM(err) *err)
{
    /* Some simple error checks */
    /* --------------------------------------------------------------------- */
    CHARM(err_check_distribution)(err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return CHARM_NAN;
    }


    if (cltmin > cltmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cltmin\" cannot be larger than \"cltmax\".");
        return CHARM_NAN;
    }


    if (n1 > pnmj->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n1\" cannot be larger than \"pnmj->nmax\".");
        return CHARM_NAN;
    }


    if (n2 > pnmj->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n2\" cannot be larger than \"pnmj->nmax\".");
        return CHARM_NAN;
    }



    if (m1 > n1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"m1\" cannot be larger than \"n1\".");
        return CHARM_NAN;
    }


    if (m2 > n2)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"m2\" cannot be larger than \"n2\".");
        return CHARM_NAN;
    }


    if (CHARM(leg_pnmj_check_ordering)(pnmj->ordering))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported value of \"pnmj->ordering\".");
        return CHARM_NAN;
    }
    /* --------------------------------------------------------------------- */






    /* Pre-computation of the trigonometric integrals */
    /* --------------------------------------------------------------------- */
    REAL *itrig = (REAL *)malloc((n1 + 1) * (n2 + 1) * sizeof(REAL));
    if (itrig == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        return CHARM_NAN;
    }


    unsigned int parity;
    unsigned int rem_m1_2 = m1 % 2;
    unsigned int rem_m2_2 = m2 % 2;
    if ((rem_m1_2 == 0) && (rem_m2_2 == 0))
    /* (Even, Even) */
        parity = 0;
    else if ((rem_m1_2 == 0) && (rem_m2_2 != 0))
    /* (Even, Odd) */
        parity = 1;
    else if ((rem_m1_2 != 0) && (rem_m2_2 == 0))
    /* (Odd, Even) */
        parity = 2;
    else /*if ((rem_m1_2 != 0) && (rem_m2_2 != 0)) */
    /* (Odd, Odd) */
        parity = 3;


    size_t i = 0;
    REAL k1d, k2d;
    REAL dclt = cltmax - cltmin;
    if (parity == 0)
    {
        for (unsigned long k1 = 0; k1 <= n1; k1++)
        {
            k1d = (REAL)k1;
            for (unsigned long k2 = 0; k2 <= n2; k2++)
            {
                k2d = (REAL)k2;
                itrig[i] = CHARM(integ_ccs)(cltmin, dclt, k1d, k2d);
                i += 1;
            }
        }
    }
    else if (parity == 1)
    {
        for (unsigned long k1 = 0; k1 <= n1; k1++)
        {
            k1d = (REAL)k1;
            for (unsigned long k2 = 0; k2 <= n2; k2++)
            {
                k2d = (REAL)k2;
                itrig[i] = CHARM(integ_css)(cltmin, dclt, k1d, k2d);
                i += 1;
            }
        }
    }
    else if (parity == 2)
    {
        for (unsigned long k1 = 0; k1 <= n1; k1++)
        {
            k1d = (REAL)k1;
            for (unsigned long k2 = 0; k2 <= n2; k2++)
            {
                k2d = (REAL)k2;
                itrig[i] = CHARM(integ_scs)(cltmin, dclt, k1d, k2d);
                i += 1;
            }
        }
    }
    else /* if (parity == 3) */
    {
        for (unsigned long k1 = 0; k1 <= n1; k1++)
        {
            k1d = (REAL)k1;
            for (unsigned long k2 = 0; k2 <= n2; k2++)
            {
                k2d = (REAL)k2;
                itrig[i] = CHARM(integ_sss)(cltmin, dclt, k1d, k2d);
                i += 1;
            }
        }
    }
    /* ..................................................................... */






    /* Computation of the integrals */
    /* --------------------------------------------------------------------- */
    /* Useful substitution */
    size_t k1_n2p1;
    unsigned long j1, j2;


    /* Initialize the value of the integral */
    REAL ip     = PREC(0.0);
    REAL ip_tmp = PREC(0.0);


    /* Loop over the Fourier coefficients of the first Legendre  function
     * inside the integral */
    for (unsigned long k1 = 0; k1 <= n1; k1++)
    {
        if (((n1 - k1) % 2) != 0)
            /* Continue, because the respective Fourier coefficient is zero by
             * definition */
            continue;


        k1_n2p1 = k1 * (n2 + 1);


        /* Loop over the Fourier coefficients of the second Legendre function
         * inside the integral */
        ip_tmp = PREC(0.0);
        for (unsigned long k2 = 0; k2 <= n2; k2++)
        {
            if (((n2 - k2) % 2) != 0)
                /* Continue, because the respective Fourier coefficient is zero
                 * by definition */
                continue;


            j2 = CHARM(leg_pnmj_k2j)(k2);
            if (pnmj->ordering == CHARM_LEG_PMNJ)
                ip_tmp += pnmj->pnmj[m2][n2 - m2][j2] * itrig[k1_n2p1 + k2];
            else if (pnmj->ordering == CHARM_LEG_PMJN)
                ip_tmp += pnmj->pnmj[m2][j2][n2 - CHARM_MAX(m2, 2 * j2)] *
                          itrig[k1_n2p1 + k2];
        }


        j1 = CHARM(leg_pnmj_k2j)(k1);
        if (pnmj->ordering == CHARM_LEG_PMNJ)
            ip += pnmj->pnmj[m1][n1 - m1][j1] * ip_tmp;
        else if (pnmj->ordering == CHARM_LEG_PMJN)
            ip += pnmj->pnmj[m1][j1][n1 - CHARM_MAX(m1, 2 * j1)] * ip_tmp;
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    free(itrig);


    return ip;
    /* --------------------------------------------------------------------- */
}
