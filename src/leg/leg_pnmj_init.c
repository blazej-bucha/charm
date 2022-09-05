/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






CHARM(pnmj) *CHARM(leg_pnmj_init)(unsigned long nmax, int nmj_order)
{
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


    /* Save the "nmj_order" value to the "CHARM(pnmj)" struct */
    /* --------------------------------------------------------------------- */
    pnmj->nmj_order = nmj_order;
    /* --------------------------------------------------------------------- */


    /* Get the size of the array to store the Fourier coefficients of Legendre
     * functions and save it to the "CHARM(pnmj)" struct */
    /* --------------------------------------------------------------------- */
    size_t npnmj = CHARM(leg_pnmj_length)(nmax);
    pnmj->npnmj = npnmj;
    /* --------------------------------------------------------------------- */


    /* Allocate the "pnmj->pnmj" array to store the Fourier coefficients of
     * Legendre functions */
    /* --------------------------------------------------------------------- */
    pnmj->pnmj = (REAL ***)malloc((nmax + 1) * sizeof(REAL **));
    if (pnmj->pnmj == NULL)
    {
        /* Memory allocation failed, so we have to deallocate all the memory
         * that has been allocated so far before we escape this function. */
        free(pnmj->pnmj); free(pnmj);
        pnmj = NULL;


        return pnmj;
    }


    if (nmj_order == CHARM_LEG_PNMJ_ORDER_MNJ)
    {
        for (unsigned long m = 0; m <= nmax; m++)
        {
            pnmj->pnmj[m] = (REAL **)malloc((nmax + 1) * sizeof(REAL *));
            if (pnmj->pnmj[m] == NULL)
            {
                /* Memory allocation failed, so we have to deallocate all the
                 * memory that has been allocated so far before we escape this
                 * function. */
                for (unsigned long m_tmp = 0; m_tmp < m; m_tmp++)
                    free(pnmj->pnmj[m]);
                free(pnmj->pnmj);
                free(pnmj);
                pnmj = NULL;


                return pnmj;
            }
        }
    }
    else if (nmj_order == CHARM_LEG_PNMJ_ORDER_MJN)
    {
        for (unsigned long m = 0; m <= nmax; m++)
        {
            pnmj->pnmj[m] = (REAL **)malloc(((nmax / 2) + 1) * sizeof(REAL *));
            if (pnmj->pnmj[m] == NULL)
            {
                /* Memory allocation failed, so we have to deallocate all the
                 * memory that has been allocated so far before we escape this
                 * function. */
                for (unsigned long m_tmp = 0; m_tmp < m; m_tmp++)
                    free(pnmj->pnmj[m]);
                free(pnmj->pnmj);
                free(pnmj);
                pnmj = NULL;


                return pnmj;
            }
        }
    }


    pnmj->pnmj[0][0] = (REAL *)calloc(npnmj, sizeof(REAL));
    if (pnmj->pnmj[0][0] == NULL)
    {
        /* Memory allocation failed, so we have to deallocate all the memory
         * that has been allocated so far before we escape this function. */
        for (unsigned long m = 0; m <= nmax; m++)
            free(pnmj->pnmj[m]);
        free(pnmj->pnmj);
        free(pnmj);
        pnmj = NULL;


        return pnmj;
    }
    /* --------------------------------------------------------------------- */


    if (nmj_order == CHARM_LEG_PNMJ_ORDER_MNJ)
    {
        /* Now set the pointers "pnmj->pnmj[m][n]" to point to the right
         * elements of the numerical array "pnmj->pnmj[0][0]" */
        /* ----------------------------------------------------------------- */
        unsigned long nj = 0;
        for (unsigned long m = 0; m <= nmax; m++)
        {
            for (unsigned long n = 0; n <= nmax; n++)
            {
                if (n < m)
                {
                    pnmj->pnmj[m][n] = NULL;
                    continue;
                }


                pnmj->pnmj[m][n] = pnmj->pnmj[0][0] + nj;
                nj += n / 2 + 1;
            }
        }
        /* ----------------------------------------------------------------- */
    }
    else if (nmj_order == CHARM_LEG_PNMJ_ORDER_MJN)
    {
        /* Now set the pointers "pnmj->pnmj[m][j]" to point to the right
         * elements of the numerical array "pnmj->pnmj[0][0]" */
        /* ----------------------------------------------------------------- */
        unsigned long jn, idx;
        jn = idx = 0;
        for (unsigned long m = 0; m <= nmax; m++)
            for (unsigned long j = 0; j <= (nmax / 2); j++)
            {
                jn                = idx - CHARM_MAX(m, 2 * j);
                pnmj->pnmj[m][j]  = pnmj->pnmj[0][0] + jn;
                idx              += nmax - CHARM_MAX(m, 2 * j) + 1;
            }
        /* ----------------------------------------------------------------- */
    }


    return pnmj;
}
