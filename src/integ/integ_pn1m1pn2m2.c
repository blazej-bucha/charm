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
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(integ_pn1m1pn2m2)(REAL cltmin, REAL cltmax,
                             unsigned long n1, unsigned long m1,
                             unsigned long n2, unsigned long m2,
                             const CHARM(pnmj) *pnmj, CHARM(err) *err)
{
    /* Some simple error checks */
    /* --------------------------------------------------------------------- */
    if (cltmin > cltmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cltmin\" cannot be larger than \"cltmax\".");


        return (ADDP(0.0) / ADDP(0.0));
    }


    if (n1 > pnmj->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n1\" cannot be larger than \"pnmj->nmax\".");


        return (ADDP(0.0) / ADDP(0.0));
    }


    if (n2 > pnmj->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n2\" cannot be larger than \"pnmj->nmax\".");


        return (ADDP(0.0) / ADDP(0.0));
    }



    if (m1 > n1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"m1\" cannot be larger than \"n1\".");


        return (ADDP(0.0) / ADDP(0.0));
    }


    if (m2 > n2)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"m2\" cannot be larger than \"n2\".");


        return (ADDP(0.0) / ADDP(0.0));
    }
    /* --------------------------------------------------------------------- */






    /* Pre-computation of the trigonometric integrals */
    /* --------------------------------------------------------------------- */
    REAL *itrig = (REAL *)malloc((n1 + 1) * (n2 + 1) * sizeof(REAL));
    if (itrig == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        return (ADDP(0.0) / ADDP(0.0));
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
    /* Get the pointers to Fourier coefficients "pnmj->pnmj" of degrees "n1"
     * and "n2" and orders "m1" and "m2" */
    REAL *pnmj_m1n1 = pnmj->pnmj[m1][n1 - m1];
    REAL *pnmj_m2n2 = pnmj->pnmj[m2][n2 - m2];


    /* Useful substitution */
    size_t k1_n2p1;


    /* Initialize the value of the integral */
    REAL ip     = ADDP(0.0);
    REAL ip_tmp = ADDP(0.0);


    /* Loop over the Fourier coefficients of the first Legendre  function
     * inside the integral */
    for (unsigned long k1 = 0; k1 <= n1; k1++)
    {
        if (((n1 - k1) % 2) != 0)
            /* Continue, because the respective Fourier  coefficient is zero by
             * definition */
            continue;


        k1_n2p1 = k1 * (n2 + 1);


        /* Loop over the Fourier coefficients of the second Legendre function
         * inside the integral */
        ip_tmp = ADDP(0.0);
        for (unsigned long k2 = 0; k2 <= n2; k2++)
        {
            if (((n2 - k2) % 2) != 0)
                /* Continue, because the respective Fourier coefficient is zero
                 * by definition */
                continue;


            ip_tmp += pnmj_m2n2[CHARM(leg_pnmj_k2j)(k2)] * itrig[k1_n2p1 + k2];
        }


        ip += pnmj_m1n1[CHARM(leg_pnmj_k2j)(k1)] * ip_tmp;
    }
    /* --------------------------------------------------------------------- */






    return ip;
}
