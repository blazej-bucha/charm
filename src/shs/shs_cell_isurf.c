/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
#include "shs_cell_isurf_coeffs.h"
#include "shs_cell_isurf_lr.h"
#include "../integ/integ_ccs.h"
#include "../integ/integ_css.h"
#include "../integ/integ_scs.h"
#include "../integ/integ_sss.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_isurf)(const CHARM(crd) *cell,
                           const CHARM(shc) *shcs1, unsigned long nmax1,
                           const CHARM(shc) *shcs2, unsigned long nmax2,
                           unsigned long nmax3, unsigned long nmax4,
                           REAL *f, CHARM(err) *err)
{
    /* Some error checks */
    /* --------------------------------------------------------------------- */
    if (cell->type != CHARM_CRD_CELLS_GRID)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->type\" must be set to "
                       "\"CHARM_CRD_CELLS_GRID\".");


        return;
    }


    if (nmax1 > shcs1->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmax1\" cannot be larger than \"shcs1->nmax\".");


        return;
    }


    if (nmax2 > shcs2->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmax2\" cannot be larger than \"shcs2->nmax\".");


        return;
    }


    if (nmax3 > nmax4)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nmax3\" cannot be larger than \"nmax4\".");


        return;
    }


    if (!CHARM(misc_is_nearly_equal)(shcs2->mu, ADDP(1.0),
                                     CHARM(glob_threshold)) ||
        !CHARM(misc_is_nearly_equal)(shcs2->r,  ADDP(1.0),
                                     CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"shcs2->mu\" and \"shcs2->r\" have to be equal to "
                       "\"1.0\".");


        return;
    }


    /* Check whether "cell->lon" is a linearly increasing array of cells. */
    size_t cell_nlon = cell->nlon;
    size_t cell_nlat = cell->nlat;
    int err_tmp = CHARM(misc_arr_chck_lin_incr)(cell->lon, 2 * cell_nlon,
                                                0, 2, CHARM(glob_threshold2),
                                                err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    if (err_tmp != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lon\" is not a linearly increasing array "
                       "of cells within the \"threshold2\".");


        return;
    }


    err_tmp = CHARM(misc_arr_chck_lin_incr)(cell->lon, 2 * cell_nlon,
                                            1, 2, CHARM(glob_threshold2), err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    if (err_tmp != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lon\" is not a linearly increasing array "
                       "of cells within the \"threshold2\".");


        return;
    }


    /* At this point, we know that the steps between cells in "cell->lon" are
     * constant.  Here, we check if the steps between the minimum longitudes
     * and the maximum longitudes are equal and, if true, store this value in
     * a separate variable (will be necessary later for the PSLR algorithm) */
    REAL dlon;
    if (cell_nlon > 1)
    {
        if (CHARM(misc_is_nearly_equal)(cell->lon[2] - cell->lon[0],
                                        cell->lon[3] - cell->lon[1],
                                        CHARM(glob_threshold)) == 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The difference \"cell->lon[2] - cell->lon[0]\" "
                           "has to be equal to "
                           "\"cell->lon[3] - cell->lon[1]\".");


            return;
        }
        dlon = cell->lon[2] - cell->lon[0];
    }
    else
        dlon = ADDP(0.0);
    /* --------------------------------------------------------------------- */







    /* --------------------------------------------------------------------- */
    int FAILURE_glob = 0;
    REAL *cnm1cnm3 = NULL;
    REAL *cnm1snm3 = NULL;
    REAL *snm1cnm3 = NULL;
    REAL *snm1snm3 = NULL;
    REAL *DELTAlon = NULL;
    /* --------------------------------------------------------------------- */






    /* Computation of the coefficients related to the potential on the
     * irregular surface */
    /* --------------------------------------------------------------------- */
    /* Allocation of the coefficients */
    size_t size = ((nmax3 + 1) * (nmax3 + 1) * (nmax1 + 1) * (nmax1 + 1));
    cnm1cnm3 = (REAL *)calloc(size, sizeof(REAL));
    if (cnm1cnm3 == NULL)
        goto FAILURE;
    cnm1snm3 = (REAL *)calloc(size, sizeof(REAL));
    if (cnm1snm3 == NULL)
        goto FAILURE;
    snm1cnm3 = (REAL *)calloc(size, sizeof(REAL));
    if (snm1cnm3 == NULL)
        goto FAILURE;
    snm1snm3 = (REAL *)calloc(size, sizeof(REAL));
    if (snm1snm3 == NULL)
        goto FAILURE;


    CHARM(shs_cell_isurf_coeffs)(shcs1, nmax1, shcs2,
                                 nmax2, nmax3, nmax4,
                                 cnm1cnm3, cnm1snm3,
                                 snm1cnm3, snm1snm3, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto FAILURE;
    }
    /* --------------------------------------------------------------------- */






    /* Synthesis of the mean values on the irregular surface */
    /* --------------------------------------------------------------------- */
    REAL lon0 = cell->lon[0];
    DELTAlon = (REAL *)malloc(cell_nlon * sizeof(REAL));
    if (DELTAlon == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    for (size_t j = 0; j < cell_nlon; j++)
        DELTAlon[j] = cell->lon[2 * j + 1] - cell->lon[2 * j];


    REAL mur = shcs1->mu / shcs1->r;
    size = (nmax1 + 1) * (nmax3 + 1);


    /* Loop over the latitude cells */
#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(cell, DELTAlon, lon0, dlon, cell_nlat, cell_nlon) \
shared(nmax1, nmax3, f, cnm1cnm3, cnm1snm3, snm1cnm3, snm1snm3) \
shared(FAILURE_glob, mur, size, err)
#endif
    {
    /* ..................................................................... */
    /* An indicator for failed memory initializations on each thread,
     * a private variable. */
    int FAILURE_priv = 0;


    REAL      *fi = NULL;
    REAL *ipt_ccs = NULL;
    REAL *ipt_css = NULL;
    REAL *ipt_scs = NULL;
    REAL *ipt_sss = NULL;


    /* Synthesized mean values for the "i"th grid row */
    fi  = (REAL *)malloc(cell_nlon * sizeof(REAL));
    if (fi == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }


    /* Arrays to store pre-computed trigonometric integrals related to
     * integrals of associated Legendre functions */
    ipt_ccs = (REAL *)malloc(size * sizeof(REAL));
    if (ipt_ccs == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }
    ipt_css = (REAL *)malloc(size * sizeof(REAL));
    if (ipt_css == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }
    ipt_scs = (REAL *)malloc(size * sizeof(REAL));
    if (ipt_scs == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }
    ipt_sss = (REAL *)malloc(size * sizeof(REAL));
    if (ipt_sss == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }


FAILURE_1_parallel:
#if CHARM_PARALLEL
#pragma omp critical
#endif
    {
        /* At this point, "FAILURE_priv" is equal to 1 if any allocation
         * failed on a particular thread.  Now, we can add "FAILURE_priv"
         * from all the threads to get the total number threads on which an
         * allocation failure occurred (the "FAILURE_glob" variable).  Note
         * that this code block is executed with one thread at a time only, so
         * "FAILURE_glob" can safely be overwritten. */
        FAILURE_glob += FAILURE_priv;
    }


    /* Now we have to wait until all the threads get here. */
#if CHARM_PARALLEL
#pragma omp barrier
#endif
    /* OK, now let's check on each thread whether there is at least one failed
     * memory allocation among the threads. */
    if (FAILURE_glob > 0)
    {
        /* Ooops, there was indeed a memory allocation failure.  So let the
         * master thread write down the error to the "err" variable. */
#if CHARM_PARALLEL
#pragma omp master
#endif
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);
        }


        /* OK, and now all threads go to the "FAILURE_2_parallel" label to
         * deallocate all the memory that might be allocated before the
         * allocation failure. */
        goto FAILURE_2_parallel;
    }
    /* ..................................................................... */


    /* The "i"th minimum and maximum co-latitudes of computation cells derived
     * from the latitudes */
    REAL clt1, clt2;


    /* Lumped coefficients */
    REAL lc00, lc01, lc10, lc11;


    /* Useful substitutions */
    unsigned int parity;
    unsigned int rem_m1_2, rem_m3_2;
    REAL dsigma, dclt;
    size_t row_idx;


    size_t idx, idx2;


#if CHARM_PARALLEL
#pragma omp for
#endif
    for (size_t i = 0; i < cell_nlat; i++)
    {
        /* Transformation of latitudes into co-latitudes */
        clt2 = PI_2 - cell->lat[2 * i];
        clt1 = PI_2 - cell->lat[2 * i + 1];


        /* Reset "fi" to zeros */
        memset(fi, 0, cell_nlon * sizeof(REAL));


        /* Reset "idx" to zero */
        idx = 0;


        /* Pre-computation of the trigonometric integrals */
        /* ----------------------------------------------------------------- */
        idx2 = 0;
        dclt = clt2 - clt1;


        {
        REAL k1d, k3d;
        for (unsigned long k1 = 0; k1 <= nmax1; k1++)
        {
            k1d = (REAL)k1;
            for (unsigned long k3 = 0; k3 <= nmax3; k3++)
            {
                k3d = (REAL)k3;


                ipt_ccs[idx2] = CHARM(integ_ccs)(clt1, dclt, k1d, k3d);
                ipt_css[idx2] = CHARM(integ_css)(clt1, dclt, k1d, k3d);
                ipt_scs[idx2] = CHARM(integ_scs)(clt1, dclt, k1d, k3d);
                ipt_sss[idx2] = CHARM(integ_sss)(clt1, dclt, k1d, k3d);


                idx2 += 1;
            }
        }
        }
        /* ----------------------------------------------------------------- */






        /* Loop over the harmonic orders of the first spherical harmonic
         * function */
        for (unsigned long m1 = 0; m1 <= nmax1; m1++)
        {
            rem_m1_2 = m1 % 2;


        /* Loop over the harmonic orders of the second spherical harmonic
         * function */
        for (unsigned long m3 = 0; m3 <= nmax3; m3++)
        {
            rem_m3_2 = m3 % 2;


            if ((rem_m1_2 == 0) && (rem_m3_2 == 0))
            /* (Even, Even) */
                parity = 0;
            else if ((rem_m1_2 == 0) && (rem_m3_2 != 0))
            /* (Even, Odd) */
                parity = 1;
            else if ((rem_m1_2 != 0) && (rem_m3_2 == 0))
            /* (Odd, Even) */
                parity = 2;
            else /*if ((rem_m1_2 != 0) && (rem_m3_2 != 0)) */
            /* (Odd, Odd) */
                parity = 3;


            /* Variables to store the lumped coefficients */
            lc00 = lc01 = lc10 = lc11 = ADDP(0.0);


            /* Computation of the lumped coefficients */
            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
            idx2 = 0;


            if (parity == 0)
            /* (Even, Even) */
            {
                for (unsigned long k1 = 0; k1 <= nmax1; k1++)
                {
                for (unsigned long k3 = 0; k3 <= nmax3; k3++)
                {
                    lc00 += ipt_ccs[idx2] * cnm1cnm3[idx];
                    lc01 += ipt_ccs[idx2] * cnm1snm3[idx];
                    lc10 += ipt_ccs[idx2] * snm1cnm3[idx];
                    lc11 += ipt_ccs[idx2] * snm1snm3[idx];


                    idx  += 1;
                    idx2 += 1;
                } /* End of the loop over "k1" */
                } /* End of the loop over "k3" */
            }
            else if (parity == 1)
            /* (Even, Odd) */
            {
                for (unsigned long k1 = 0; k1 <= nmax1; k1++)
                {
                for (unsigned long k3 = 0; k3 <= nmax3; k3++)
                {
                    lc00 += ipt_css[idx2] * cnm1cnm3[idx];
                    lc01 += ipt_css[idx2] * cnm1snm3[idx];
                    lc10 += ipt_css[idx2] * snm1cnm3[idx];
                    lc11 += ipt_css[idx2] * snm1snm3[idx];


                    idx  += 1;
                    idx2 += 1;
                } /* End of the loop over "k1" */
                } /* End of the loop over "k3" */
            }
            else if (parity == 2)
            /* (Odd, Even) */
            {
                for (unsigned long k1 = 0; k1 <= nmax1; k1++)
                {
                for (unsigned long k3 = 0; k3 <= nmax3; k3++)
                {
                    lc00 += ipt_scs[idx2] * cnm1cnm3[idx];
                    lc01 += ipt_scs[idx2] * cnm1snm3[idx];
                    lc10 += ipt_scs[idx2] * snm1cnm3[idx];
                    lc11 += ipt_scs[idx2] * snm1snm3[idx];


                    idx  += 1;
                    idx2 += 1;
                } /* End of the loop over "k1" */
                } /* End of the loop over "k3" */
            }
            else /*if (parity == 3) */
            /* (Odd, Odd) */
            {
                for (unsigned long k1 = 0; k1 <= nmax1; k1++)
                {
                for (unsigned long k3 = 0; k3 <= nmax3; k3++)
                {
                    lc00 += ipt_sss[idx2] * cnm1cnm3[idx];
                    lc01 += ipt_sss[idx2] * cnm1snm3[idx];
                    lc10 += ipt_sss[idx2] * snm1cnm3[idx];
                    lc11 += ipt_sss[idx2] * snm1snm3[idx];


                    idx  += 1;
                    idx2 += 1;
                } /* End of the loop over "k1" */
                } /* End of the loop over "k3" */
            }
            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


            /* Synthesis for the "i"th latitude parallel */
            CHARM(shs_cell_isurf_lr)(lon0, dlon, cell_nlon,
                                     lc00, lc01, lc10, lc11, m1, m3,
                                     fi);
        }  /* End of the loop over "m1" */
        }  /* End of the loop over "m3" */


        dclt = COS(clt1) - COS(clt2);
        row_idx = i * cell_nlon;
        for (size_t j = 0; j < cell_nlon; j++)
        {
            /* Compute the area of the cells on the unit sphere */
            dsigma = dclt * DELTAlon[j];


            /* Final synthesis */
            f[row_idx + j] = mur / dsigma * fi[j];
        }


    }  /* End of the loop over "i" */


    /* Free the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE_2_parallel:
    free(fi);
    free(ipt_ccs);
    free(ipt_css);
    free(ipt_scs);
    free(ipt_sss);
    /* --------------------------------------------------------------------- */


    }
    /* --------------------------------------------------------------------- */






    /* Free the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if (FAILURE_glob != 0)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(DELTAlon);
    free(cnm1cnm3); free(cnm1snm3); free(snm1cnm3); free(snm1snm3);
    /* --------------------------------------------------------------------- */





    return;
}
