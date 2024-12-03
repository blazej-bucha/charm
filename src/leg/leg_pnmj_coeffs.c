/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include "../prec.h"
#include "leg_pnmj_leq4.h"
#include "leg_pnmj_dpeven.h"
#include "leg_pnmj_gpeven.h"
#include "leg_pnmj_dpodd.h"
#include "leg_pnmj_gpodd.h"
#include "leg_pnmj_check_ordering.h"
#include "../xnum/xnum_x2f.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_check_distribution.h"
/* ------------------------------------------------------------------------- */






void CHARM(leg_pnmj_coeffs)(CHARM(pnmj) *pnmj,
                            unsigned long nmax,
                            CHARM(err) *err)
{
    /* Error checks */
    /* --------------------------------------------------------------------- */
    CHARM(err_check_distribution)(err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    if (nmax > pnmj->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmax\" cannot be larger than \"pnmj->nmax\".");
        return;
    }


    if (CHARM(leg_pnmj_check_ordering)(pnmj->ordering))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported value of \"pnmj->ordering\".");
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Initializations */
    /* --------------------------------------------------------------------- */
    unsigned long nmaxp1 = nmax + 1;


    REAL *xp    = NULL;
    REAL *xpold = NULL;
    REAL *xp0   = NULL;
    REAL *xp1   = NULL;
    REAL *xp2   = NULL;


    int *ip    = NULL;
    int *ipold = NULL;
    int *ip0   = NULL;
    int *ip1   = NULL;
    int *ip2   = NULL;


    xp    = (REAL *)calloc(nmaxp1, sizeof(REAL));
    if (xp == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    xpold = (REAL *)calloc(nmaxp1, sizeof(REAL));
    if (xpold == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    xp0   = (REAL *)calloc(nmaxp1, sizeof(REAL));
    if (xp0 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    xp1   = (REAL *)calloc(nmaxp1, sizeof(REAL));
    if (xp1 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    xp2   = (REAL *)calloc(nmaxp1, sizeof(REAL));
    if (xp2 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }


    ip    = (int *)calloc(nmaxp1, sizeof(int));
    if (ip == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    ipold = (int *)calloc(nmaxp1, sizeof(int));
    if (ipold == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    ip0   = (int *)calloc(nmaxp1, sizeof(int));
    if (ip0 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    ip1   = (int *)calloc(nmaxp1, sizeof(int));
    if (ip1 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    ip2   = (int *)calloc(nmaxp1, sizeof(int));
    if (ip2 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }


    /* Before doing any actual computations, reset all Fourier coefficients to
     * zero.  This ensures that if "nmax < pnmj->nmax", all Fourier
     * coefficients in "pnmj->pnmj" beyond "nmax" will be set to zero. */
    memset(pnmj->pnmj[0][0], 0, pnmj->npnmj * sizeof(REAL));
    /* --------------------------------------------------------------------- */






    /* Even harmonic degrees up to "n == 4" */
    /* --------------------------------------------------------------------- */
    REAL tmp1, tmp2;
    unsigned long jmax, nmaxtmp = CHARM_MIN(nmax, 4);
    for (unsigned long n = 0; n <= nmaxtmp; n += 2)
    {
        jmax = n / 2;


        for (unsigned long m = 0; m <= n; m++)
        {
            CHARM(leg_pnmj_leq4)(n, m, xp, err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto FAILURE;
            }


            for (unsigned long j = 0; j <= jmax; j++)
            {
                ip[j] = 0;


                tmp1 = CHARM(xnum_x2f)(xp[j], ip[j]);
                if (pnmj->ordering == CHARM_LEG_PMNJ)
                    pnmj->pnmj[m][n - m][j] = tmp1;
                else if (pnmj->ordering == CHARM_LEG_PMJN)
                    pnmj->pnmj[m][j][n - CHARM_MAX(m, 2 * j)] = tmp1;
            }
        }

    }
    /* --------------------------------------------------------------------- */






    /* Even harmonic degrees from "n = 6" up to "n == nmax" */
    /* --------------------------------------------------------------------- */
    for (unsigned long n = 6; n <= nmax; n += 2)
    {
        for (unsigned long j = 0; j <= jmax; j++)
        {
            xpold[j] = xp[j];
            ipold[j] = ip[j];
        }


        CHARM(leg_pnmj_dpeven)(n, xpold, xp, xp1, ipold, ip, ip1, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto FAILURE;
        }


        jmax = n / 2;
        for (unsigned long j = 0; j <= jmax; j++)
        {
            tmp1 = CHARM(xnum_x2f)(xp[j],  ip[j]);
            tmp2 = CHARM(xnum_x2f)(xp1[j], ip1[j]);
            if (pnmj->ordering == CHARM_LEG_PMNJ)
            {
                pnmj->pnmj[n][0][j]     = tmp1;
                pnmj->pnmj[n - 1][1][j] = tmp2;
            }
            else if (pnmj->ordering == CHARM_LEG_PMJN)
            {
                pnmj->pnmj[n][j][0]                               = tmp1;
                pnmj->pnmj[n - 1][j][n - CHARM_MAX(2 * j, n - 1)] = tmp2;
            }


            xp2[j] = xp[j];
            ip2[j] = ip[j];
        }


        /* The loop can alternatively be designed as
         *
         *      for (unsigned long m = n - 2; m != (unsigned long)-1; m--)
         *
         * but the current form seems to be more robust. */
        for (unsigned long m = n - 1; m-- > 0;)
        {
            CHARM(leg_pnmj_gpeven)(jmax, n, m, xp2, xp1, xp0, ip2, ip1,
                                   ip0, err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto FAILURE;
            }


            for (unsigned long j = 0; j <= jmax; j++)
            {
                tmp1 = CHARM(xnum_x2f)(xp0[j], ip0[j]);
                if (pnmj->ordering == CHARM_LEG_PMNJ)
                    pnmj->pnmj[m][n - m][j] = tmp1;
                else if (pnmj->ordering == CHARM_LEG_PMJN)
                    pnmj->pnmj[m][j][n - CHARM_MAX(m, 2 * j)] = tmp1;


                xp2[j] = xp1[j];
                ip2[j] = ip1[j];


                xp1[j] = xp0[j];
                ip1[j] = ip0[j];
            }
        }
    }
    /* --------------------------------------------------------------------- */






    /* Odd harmonic degrees up to "n == 3" */
    /* --------------------------------------------------------------------- */
    nmaxtmp = CHARM_MIN(nmax, 3);
    for (unsigned long n = 1; n <= nmaxtmp; n += 2)
    {
        jmax = (n - 1) / 2;


        for (unsigned long m = 0; m <= n; m++)
        {
            CHARM(leg_pnmj_leq4)(n, m, xp, err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto FAILURE;
            }


            for (unsigned long j = 0; j <= jmax; j++)
            {
                ip[j] = 0;
                tmp1 = CHARM(xnum_x2f)(xp[j], ip[j]);
                if (pnmj->ordering == CHARM_LEG_PMNJ)
                    pnmj->pnmj[m][n - m][j] = tmp1;
                else if (pnmj->ordering == CHARM_LEG_PMJN)
                    pnmj->pnmj[m][j][n - CHARM_MAX(m, 2 * j)] = tmp1;
            }
        }
    }
    /* --------------------------------------------------------------------- */






    /* Odd harmonic degrees from "n == 5" up to "n == nmax" */
    /* --------------------------------------------------------------------- */
    for (unsigned long n = 5; n <= nmax; n += 2)
    {
        for (unsigned long j = 0; j <= jmax; j++)
        {
            xpold[j] = xp[j];
            ipold[j] = ip[j];
        }


        CHARM(leg_pnmj_dpodd)(n, xpold, xp, xp1, ipold, ip, ip1, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto FAILURE;
        }


        jmax = (n - 1) / 2;
        for (unsigned long j = 0; j <= jmax; j++)
        {
            tmp1 = CHARM(xnum_x2f)(xp[j],  ip[j]);
            tmp2 = CHARM(xnum_x2f)(xp1[j], ip1[j]);
            if (pnmj->ordering == CHARM_LEG_PMNJ)
            {
                pnmj->pnmj[n][0][j]     = tmp1;
                pnmj->pnmj[n - 1][1][j] = tmp2;
            }
            else if (pnmj->ordering == CHARM_LEG_PMJN)
            {
                pnmj->pnmj[n][j][0]                               = tmp1;
                pnmj->pnmj[n - 1][j][n - CHARM_MAX(2 * j, n - 1)] = tmp2;
            }


            xp2[j] = xp[j];
            ip2[j] = ip[j];
        }


        /* The loop can alternatively be designed as
         *
         *      for (unsigned long m = n - 2; m != (unsigned long)-1; m--)
         *
         * but the current form seems to be more robust. */
        for (unsigned long m = n - 1; m-- > 0;)
        {
            CHARM(leg_pnmj_gpodd)(jmax, n, m, xp2, xp1, xp0, ip2, ip1, ip0,
                                    err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto FAILURE;
            }


            for (unsigned long j = 0; j <= jmax; j++)
            {
                tmp1 = CHARM(xnum_x2f)(xp0[j], ip0[j]);
                if (pnmj->ordering == CHARM_LEG_PMNJ)
                    pnmj->pnmj[m][n - m][j] = tmp1;
                else if (pnmj->ordering == CHARM_LEG_PMJN)
                    pnmj->pnmj[m][j][n - CHARM_MAX(m, 2 * j)] = tmp1;


                xp2[j] = xp1[j];
                ip2[j] = ip1[j];


                xp1[j] = xp0[j];
                ip1[j] = ip0[j];
            }
        }
    }
    /* --------------------------------------------------------------------- */






    /* Free the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    free(xp); free(xpold); free(xp0); free(xp1); free(xp2);
    free(ip); free(ipold); free(ip0); free(ip1); free(ip2);
    /* --------------------------------------------------------------------- */






    return;
}
