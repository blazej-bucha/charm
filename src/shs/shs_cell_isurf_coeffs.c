/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "shs_cell_isurf_offset.h"
#include "shs_cell_isurf_coeffs.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_isurf_coeffs)(const CHARM(shc) *shcs1, unsigned long nmax1,
                                  const CHARM(shc) *shcs2, unsigned long nmax2,
                                  unsigned long nmax3,
                                  unsigned long nmax4,
                                  REAL *cnm1cnm3, REAL *cnm1snm3,
                                  REAL *snm1cnm3, REAL *snm1snm3,
                                  CHARM(err) *err)
/*
 * ============================================================================
 *
 * DESCRIPTION: Prepares the coefficients of the potential on the irregular
 *              surface to be used in "CHARM(shs_cell_isurf)".
 *
 *
 * INPUTS: "shcs1", ..., "nmax4" -- All variables have the same meaning
 *                                  as in "CHARM(shs_cell_isurf)", where
 *                                   further details can be found.
 *
 *
 * OUTPUTS: "cnm1cnm3" -- A pointer to an array of dimensions "(nmax1 + 1) *
 *                        (nmax1 + 1) * (nmax3 + 1) * (nmax3 + 1)" to store the
 *                        coefficients.  The array has 4 dimensions and should
 *                        be allocated and initialized to zeros before this
 *                        functions is called as follows:
 *
 *                        double *cnm1cnm3 = (double *)calloc(((nmax1 + 1)
 *                              * (nmax1 + 1) * (nmax3 + 1) * (nmax3 + 1)),
 *                              sizeof(double));
 *
 *                        The structure of the array can be easily guessed from
 *                        how the "lc" variable is treated here and combined
 *                        with the "CHARM(shs_cell_isurf_offset)" function.
 *
 *          "cnm1snm3", "snm1cnm3", "snm1snm3" -- Same as "cnm1cnm3", but for a
 *                        different combination of "sine" and "cosine" terms.
 *
 * ============================================================================
 *
 * */
{
    /* --------------------------------------------------------------------- */
    CHARM(point) *glg     = NULL;
    CHARM(shc) *shcs3     = NULL;
    CHARM(pnmj) *pnmj     = NULL;
    CHARM(pnmj) *cnm1pnmj = NULL;
    CHARM(pnmj) *snm1pnmj = NULL;
    REAL *r               = NULL;
    REAL *r_pow           = NULL;
    REAL *cnm3pnmj_sum    = NULL;
    REAL *snm3pnmj_sum    = NULL;
    /* --------------------------------------------------------------------- */






    /* Prepare the Gauss--Legendre grid */
    /* --------------------------------------------------------------------- */
    glg = CHARM(crd_point_gl)(nmax4, PREC(1.0));
    if (glg == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       "Failed to initialize the Gauss--Legendre grid.");
        goto FAILURE;
    }
    /* --------------------------------------------------------------------- */






    /* Synthesize the topography at the nodes of the Gauss--Legendre grid */
    /* --------------------------------------------------------------------- */
    r = (REAL *)malloc(glg->npoint * sizeof(REAL));
    if (r == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }


    CHARM(shs_point)(glg, shcs2, nmax2, r, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto FAILURE;
    }


    /* "shcs1->r / r" */
    size_t i;
#if CHARM_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i = 0; i < glg->npoint; i++)
        r[i] = shcs1->r / r[i];


    /* Array to store the "n + 1"th power of "shcs1->r / r" for "n = 0, 1, ...,
     * nmax1" */
    r_pow = (REAL *)malloc(glg->npoint * sizeof(REAL));
    if (r_pow == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
#if CHARM_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
    for (i = 0; i < glg->npoint; i++)
        r_pow[i] = PREC(1.0);
    /* --------------------------------------------------------------------- */






    /* Initialize the spherical harmonic coefficients for the "n + 1"th power
     * */
    /* --------------------------------------------------------------------- */
    shcs3 = CHARM(shc_calloc)(nmax3, PREC(1.0), PREC(1.0));
    if (shcs3 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       "Failed to initialize the \"shc\" structure.");
        goto FAILURE;
    }


    REAL *cnm3_m3, *snm3_m3;
    /* --------------------------------------------------------------------- */






    /* Pre-computation of the Fourier coefficients of Legendre functions */
    /* --------------------------------------------------------------------- */
    unsigned long nmax = CHARM_MAX(nmax1, nmax3);


    pnmj = CHARM(leg_pnmj_calloc)(nmax, CHARM_LEG_PMJN);
    if (pnmj == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       "Failed to initialize the \"pnmj\" structure.");
        goto FAILURE;
    }


    CHARM(leg_pnmj_coeffs)(pnmj, nmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto FAILURE;
    }


    REAL *pnmj_m1_j1, *pnmj_m3_j3;
    /* --------------------------------------------------------------------- */






    /* Useful substitutions */
    /* --------------------------------------------------------------------- */
    unsigned long n1_2, n1_rem_2, j1pj1, j3pj3;
    unsigned long nmax3_2 = nmax3 / 2;
    unsigned long nmax1p1 = nmax1 + 1, nmax3p1 = nmax3 + 1;
    unsigned long max_m1_j1, max_m3_j3;
    unsigned long k1, k3, k3pm1;


    /* Variables to store pre-computed products of "shcs1->c", "shcs1->s" and
     * "pnmj->pnmj". It is useful to use the same structure as for "pnmj";
     * therefore, the "CHARM(leg_pnmj_calloc)" function can be employed to
     * initialize the two variables. */
    /* ..................................................................... */
    cnm1pnmj = CHARM(leg_pnmj_calloc)(nmax1, CHARM_LEG_PMJN);
    if (cnm1pnmj == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       "Failed to initialize the \"pnmj\" structure.");
        goto FAILURE;
    }


    snm1pnmj = CHARM(leg_pnmj_calloc)(nmax1, CHARM_LEG_PMJN);
    if (snm1pnmj == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       "Failed to initialize the \"pnmj\" structure.");
        goto FAILURE;
    }
    /* ..................................................................... */


    REAL **cnm1pnmj_m1, *cnm1pnmj_m1_j1, cnm1pnmj_m1_j1_n1;
    REAL **snm1pnmj_m1, *snm1pnmj_m1_j1, snm1pnmj_m1_j1_n1;


    /* Variables to store pre-computed sums of the products "shcs3->c",
     * "shcs3->s" and "pnmj->pnmj" */
    cnm3pnmj_sum = (REAL *)malloc(nmax3p1 * (nmax3_2 + 1) * 2 *
                                    sizeof(REAL));
    if (cnm3pnmj_sum == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    snm3pnmj_sum = (REAL *)malloc(nmax3p1 * (nmax3_2 + 1) * 2 *
                                    sizeof(REAL));
    if (snm3pnmj_sum == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    REAL c_sum1, c_sum2, s_sum1, s_sum2;


    size_t idx, idxp1, idx0_lc, idx1_lc;
    /* --------------------------------------------------------------------- */






    /* The actual computation of the coefficients */
    /* --------------------------------------------------------------------- */
    unsigned long m1;
#if CHARM_OPENMP
#pragma omp parallel for default(none) \
shared(shcs1, cnm1pnmj, snm1pnmj, pnmj, nmax1) \
private(cnm1pnmj_m1, cnm1pnmj_m1_j1, snm1pnmj_m1, snm1pnmj_m1_j1) \
private(m1, pnmj_m1_j1, max_m1_j1)
#endif
    for (m1 = 0; m1 <= nmax1; m1++)
    {
        cnm1pnmj_m1 = cnm1pnmj->pnmj[m1];
        snm1pnmj_m1 = snm1pnmj->pnmj[m1];


        for (unsigned long j1 = 0; j1 <= (nmax1 / 2); j1++)
        {
            max_m1_j1 = CHARM_MAX(m1, j1 + j1);


            cnm1pnmj_m1_j1 = cnm1pnmj->pnmj[m1][j1];
            snm1pnmj_m1_j1 = snm1pnmj->pnmj[m1][j1];
            pnmj_m1_j1     = pnmj->pnmj[m1][j1];


            for (unsigned long n1 = max_m1_j1; n1 <= nmax1; n1++)
            {
                cnm1pnmj_m1_j1[n1 - max_m1_j1] = shcs1->c[m1][n1 - m1] *
                                                 pnmj_m1_j1[n1 - max_m1_j1];
                snm1pnmj_m1_j1[n1 - max_m1_j1] = shcs1->s[m1][n1 - m1] *
                                                 pnmj_m1_j1[n1 - max_m1_j1];
            }
        }
    }


    for (unsigned long n1 = 0; n1 <= nmax1; n1++)
    {
        n1_2 = n1 / 2;
        n1_rem_2 = n1 % 2;


        /* Compute the "n1 + 1"th power of the topographic surface and its
         * spherical harmonic coefficients */
        /* ................................................................. */
        size_t i;
#if CHARM_OPENMP
#pragma omp parallel for default(shared) private(i)
#endif
        for (i = 0; i < glg->npoint; i++)
            r_pow[i] *= r[i];


        CHARM(sha_point)(glg, r_pow, nmax3, shcs3, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto FAILURE;
        }
        /* ................................................................. */


        /* Pre-compute the summation over "n3" */
        /* ................................................................. */
        unsigned long m3;
#if CHARM_OPENMP
#pragma omp parallel for default(none) \
shared(nmax3, nmax3_2, pnmj, shcs3, cnm3pnmj_sum, snm3pnmj_sum) \
private(j3pj3, max_m3_j3, idx, pnmj_m3_j3, c_sum1, c_sum2, s_sum1, s_sum2) \
private(m3, cnm3_m3, snm3_m3)
#endif
        for (m3 = 0; m3 <= nmax3; m3++)
        {
            cnm3_m3 = shcs3->c[m3];
            snm3_m3 = shcs3->s[m3];


            for (unsigned long j3 = 0; j3 <= nmax3_2; j3++)
            {
                j3pj3     = j3 + j3;
                max_m3_j3 = CHARM_MAX(m3, j3pj3);
                idx       = (m3 * (nmax3_2 + 1) + j3) * 2;


                c_sum1 = c_sum2 = PREC(0.0);
                s_sum1 = s_sum2 = PREC(0.0);
                pnmj_m3_j3      = pnmj->pnmj[m3][j3];


                {
                unsigned long n3;
                for (n3 = max_m3_j3; n3 < nmax3; n3 += 2)
                {
                    c_sum1 += cnm3_m3[n3 - m3] * pnmj_m3_j3[n3 - max_m3_j3];
                    s_sum1 += snm3_m3[n3 - m3] * pnmj_m3_j3[n3 - max_m3_j3];


                    c_sum2 += cnm3_m3[n3 + 1 - m3] *
                              pnmj_m3_j3[n3 + 1 - max_m3_j3];
                    s_sum2 += snm3_m3[n3 + 1 - m3] *
                              pnmj_m3_j3[n3 + 1 - max_m3_j3];
                }


                if (n3 == nmax3)
                {
                    c_sum1 += cnm3_m3[n3 - m3] * pnmj_m3_j3[n3 - max_m3_j3];
                    s_sum1 += snm3_m3[n3 - m3] * pnmj_m3_j3[n3 - max_m3_j3];
                }
                }


                cnm3pnmj_sum[idx]     = c_sum1;
                cnm3pnmj_sum[idx + 1] = c_sum2;
                snm3pnmj_sum[idx]     = s_sum1;
                snm3pnmj_sum[idx + 1] = s_sum2;
            }
        }
        /* ................................................................. */


        /* Compute the contribution of "n1" to coefficients and add it to the
         * contribution from "n1 - 1" */
        /* ................................................................. */
        unsigned long m1;
#if CHARM_OPENMP
#pragma omp parallel for default(none) \
shared(n1, n1_2, n1_rem_2, nmax3, nmax3_2, nmax1p1, nmax3p1) \
shared(cnm1pnmj, snm1pnmj, cnm3pnmj_sum, snm3pnmj_sum) \
shared(cnm1cnm3, cnm1snm3, snm1cnm3, snm1snm3) \
private(cnm1pnmj_m1, snm1pnmj_m1, cnm1pnmj_m1_j1, snm1pnmj_m1_j1) \
private(cnm1pnmj_m1_j1_n1, snm1pnmj_m1_j1_n1) \
private(m1, k1, k3, k3pm1, j1pj1, j3pj3, max_m1_j1, max_m3_j3) \
private(idx, idxp1, idx0_lc, idx1_lc)
#endif
        for (m1 = 0; m1 <= n1; m1++)
        {
            cnm1pnmj_m1 = cnm1pnmj->pnmj[m1];
            snm1pnmj_m1 = snm1pnmj->pnmj[m1];


            for (unsigned long m3 = 0; m3 <= nmax3; m3++)
            {
                for (unsigned long j1 = 0; j1 <= n1_2; j1++)
                {
                    j1pj1 = j1 + j1;
                    max_m1_j1 = CHARM_MAX(m1, j1pj1);
                    k1 = n1_rem_2 + j1pj1;


                    cnm1pnmj_m1_j1_n1 = cnm1pnmj_m1[j1][n1 - max_m1_j1];
                    snm1pnmj_m1_j1_n1 = snm1pnmj_m1[j1][n1 - max_m1_j1];


                    for (unsigned long j3 = 0; j3 <= nmax3_2; j3++)
                    {
                        j3pj3     = j3 + j3;
                        max_m3_j3 = CHARM_MAX(m3, j3pj3);
                        k3        = (max_m3_j3 % 2) + j3pj3;


                        idx     = (m3 * (nmax3_2 + 1) + j3) * 2;
                        idx0_lc = CHARM(shs_cell_isurf_offset)(m1, m3, k1, k3,
                                                            nmax3p1, nmax1p1,
                                                            nmax3p1);


                        cnm1cnm3[idx0_lc] += cnm1pnmj_m1_j1_n1 *
                                             cnm3pnmj_sum[idx];
                        cnm1snm3[idx0_lc] += cnm1pnmj_m1_j1_n1 *
                                             snm3pnmj_sum[idx];
                        snm1cnm3[idx0_lc] += snm1pnmj_m1_j1_n1 *
                                             cnm3pnmj_sum[idx];
                        snm1snm3[idx0_lc] += snm1pnmj_m1_j1_n1 *
                                             snm3pnmj_sum[idx];


                        if (max_m3_j3 < nmax3)
                        {
                            idxp1   = idx + 1;
                            k3pm1   = (k3 % 2) ? k3 - 1 : k3 + 1;
                            idx1_lc = CHARM(shs_cell_isurf_offset)(m1, m3, k1,
                                                               k3pm1,
                                                               nmax3p1,
                                                               nmax1p1,
                                                               nmax3p1);


                            cnm1cnm3[idx1_lc] += cnm1pnmj_m1_j1_n1 *
                                                 cnm3pnmj_sum[idxp1];
                            cnm1snm3[idx1_lc] += cnm1pnmj_m1_j1_n1 *
                                                 snm3pnmj_sum[idxp1];
                            snm1cnm3[idx1_lc] += snm1pnmj_m1_j1_n1 *
                                                 cnm3pnmj_sum[idxp1];
                            snm1snm3[idx1_lc] += snm1pnmj_m1_j1_n1 *
                                                 snm3pnmj_sum[idxp1];
                        }
                    }
                }
            }
        }
        /* ................................................................. */
    }
    /* --------------------------------------------------------------------- */






    /* Free the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    CHARM(crd_point_free)(glg);
    free(r);
    free(r_pow);
    free(cnm3pnmj_sum);
    free(snm3pnmj_sum);
    CHARM(leg_pnmj_free)(pnmj);
    CHARM(leg_pnmj_free)(cnm1pnmj);
    CHARM(leg_pnmj_free)(snm1pnmj);
    CHARM(shc_free)(shcs3);
    /* --------------------------------------------------------------------- */






    return;
}
