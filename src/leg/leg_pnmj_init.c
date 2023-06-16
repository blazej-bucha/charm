/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
#include "leg_pnmj_length.h"
#include "leg_pnmj_check_ordering.h"
/* ------------------------------------------------------------------------- */






CHARM(pnmj) *CHARM(leg_pnmj_init)(unsigned long nmax,
                                  int ordering,
                                  REAL *pnmj_coeffs)
{
    /* Check inputs */
    /* --------------------------------------------------------------------- */
    if (CHARM(leg_pnmj_check_ordering)(ordering))
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Allocate memory for the "CHARM(pnmj)" data type */
    /* --------------------------------------------------------------------- */
    CHARM(pnmj) *pnmj = (CHARM(pnmj) *)malloc(sizeof(CHARM(pnmj)));
    if (pnmj == NULL)
        return pnmj;
    /* --------------------------------------------------------------------- */


    /* Save the "nmax" value to the "CHARM(pnmj)" struct */
    /* --------------------------------------------------------------------- */
    pnmj->nmax = nmax;
    /* --------------------------------------------------------------------- */


    /* Save the "ordering" value to the "CHARM(pnmj)" struct */
    /* --------------------------------------------------------------------- */
    pnmj->ordering = ordering;
    /* --------------------------------------------------------------------- */


    /* Get the size of the array to store the Fourier coefficients of Legendre
     * functions and save it to the "CHARM(pnmj)" struct */
    /* --------------------------------------------------------------------- */
    pnmj->npnmj = CHARM(leg_pnmj_length)(nmax);
    /* --------------------------------------------------------------------- */


    /* Allocate the "pnmj->pnmj" array to store the Fourier coefficients of
     * Legendre functions */
    /* --------------------------------------------------------------------- */
    pnmj->pnmj = (REAL ***)malloc((nmax + 1) * sizeof(REAL **));
    if (pnmj->pnmj == NULL)
    {
        /* Memory allocation failed, so we have to deallocate all the memory
         * that has been allocated so far before we escape this function. */
        free(pnmj->pnmj);
        free(pnmj);
        return NULL;
    }


    if (ordering == CHARM_LEG_PMNJ)
    {
        for (unsigned long m = 0; m <= nmax; m++)
        {
            pnmj->pnmj[m] = (REAL **)malloc((nmax + 1 - m) * sizeof(REAL *));
            if (pnmj->pnmj[m] == NULL)
            {
                /* Memory allocation failed, so we have to deallocate all the
                 * memory that has been allocated so far before we escape this
                 * function. */
                for (unsigned long m_tmp = 0; m_tmp < m; m_tmp++)
                    free(pnmj->pnmj[m]);
                free(pnmj->pnmj);
                free(pnmj);
                return NULL;
            }
        }
    }
    else if (ordering == CHARM_LEG_PMJN)
    {
        for (unsigned long m = 0; m <= nmax; m++)
        {
            pnmj->pnmj[m] = (REAL **)malloc(((nmax / 2) + 1) *
                                              sizeof(REAL *));
            if (pnmj->pnmj[m] == NULL)
            {
                /* Memory allocation failed, so we have to deallocate all the
                 * memory that has been allocated so far before we escape this
                 * function. */
                for (unsigned long m_tmp = 0; m_tmp < m; m_tmp++)
                    free(pnmj->pnmj[m]);
                free(pnmj->pnmj);
                free(pnmj);
                return NULL;
            }
        }
    }


    /* Assign the "pnmj_coefficients" to the "charm_pnmj" structure. */
    pnmj->pnmj[0][0] = pnmj_coeffs;
    /* --------------------------------------------------------------------- */


    if (ordering == CHARM_LEG_PMNJ)
    {
        /* Now set the pointers "pnmj->pnmj[m][n - m]" to point to the right
         * elements of the numerical array "pnmj->pnmj[0][0]" */
        /* ----------------------------------------------------------------- */
        unsigned long nj = 0;
        for (unsigned long m = 0; m <= nmax; m++)
            for (unsigned long n = m; n <= nmax; n++)
            {
                pnmj->pnmj[m][n - m] = pnmj->pnmj[0][0] + nj;
                nj += (n / 2) + 1;
            }
        /* ----------------------------------------------------------------- */
    }
    else if (ordering == CHARM_LEG_PMJN)
    {
        /* Now set the pointers "pnmj->pnmj[m][j]" to point to the right
         * elements of the numerical array "pnmj->pnmj[0][0]" */
        /* ----------------------------------------------------------------- */
        unsigned long jn = 0;
        for (unsigned long m = 0; m <= nmax; m++)
            for (unsigned long j = 0; j <= (nmax / 2); j++)
            {
                pnmj->pnmj[m][j] = pnmj->pnmj[0][0] + jn;
                jn += nmax - CHARM_MAX(2 * j, m) + 1;
            }
        /* ----------------------------------------------------------------- */
    }


    return pnmj;
}
