/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../prec.h"
#include "shs_cell_isurf_coeffs.h"
#include "shs_cell_isurf_lr.h"
#include "shs_cell_check_grd_lons.h"
#include "../integ/integ_ccs.h"
#include "../integ/integ_css.h"
#include "../integ/integ_scs.h"
#include "../integ/integ_sss.h"
#include "../crd/crd_check_cells.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../simd/simd.h"
#include "../simd/malloc_aligned.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
#include "../misc/misc_is_nearly_equal.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_isurf)(const CHARM(cell) *cell,
                           const CHARM(shc) *shcs1, unsigned long nmax1,
                           const CHARM(shc) *shcs2, unsigned long nmax2,
                           unsigned long nmax3, unsigned long nmax4,
                           REAL *f, CHARM(err) *err)
{
    /* Some error checks */
    /* --------------------------------------------------------------------- */
    if (cell->type != CHARM_CRD_CELL_GRID)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->type\" must be set to "
                       "\"CHARM_CRD_CELL_GRID\".");
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


    if (!CHARM(misc_is_nearly_equal)(shcs2->mu, PREC(1.0),
                                     CHARM(glob_threshold)) ||
        !CHARM(misc_is_nearly_equal)(shcs2->r,  PREC(1.0),
                                     CHARM(glob_threshold)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"shcs2->mu\" and \"shcs2->r\" have to be equal to "
                       "\"1.0\".");
        return;
    }


    /* Check whether "cell->lonmin" and "cell->lonmax" are linearly increasing
     * arrays. */
    size_t cell_nlon = cell->nlon;
    size_t cell_nlat = cell->nlat;
    REAL dlon;
    CHARM(shs_cell_check_grd_lons)(cell, &dlon, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    /* Check cell boundaries */
    CHARM(crd_check_cells)(cell, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
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
    size_t size = (nmax3 + 1) * (nmax3 + 1) * (nmax1 + 1) * (nmax1 + 1);
    cnm1cnm3 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, size,
                                             sizeof(REAL));
    if (cnm1cnm3 == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    cnm1snm3 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, size,
                                             sizeof(REAL));
    if (cnm1snm3 == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    snm1cnm3 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, size,
                                             sizeof(REAL));
    if (snm1cnm3 == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    snm1snm3 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, size,
                                             sizeof(REAL));
    if (snm1snm3 == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }


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






    /* Synthesis of the mean values on the irregular surface.  No polar
     * optimization is used here. */
    /* --------------------------------------------------------------------- */
    REAL lon0 = cell->lonmin[0];
    DELTAlon = (REAL *)malloc(cell_nlon * sizeof(REAL));
    if (DELTAlon == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    for (size_t j = 0; j < cell_nlon; j++)
        DELTAlon[j] = cell->lonmax[j] - cell->lonmin[j];


    REAL mur = shcs1->mu / shcs1->r;
    size = (nmax1 + 1) * (nmax3 + 1) * SIMD_SIZE;


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
    REAL   *clt1v = NULL;
    REAL   *clt2v = NULL;
    REAL   *dcltv = NULL;


    /* Synthesized mean values for the "i"th grid row */
    fi = (REAL *)CHARM(malloc_aligned)(SIMD_MEMALIGN, cell_nlon * SIMD_SIZE *
                                       sizeof(REAL));
    if (fi == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }


    /* Arrays to store pre-computed trigonometric integrals related to
     * integrals of associated Legendre functions */
    ipt_ccs = (REAL *)CHARM(malloc_aligned)(SIMD_MEMALIGN,
                                            size * sizeof(REAL));
    if (ipt_ccs == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }
    ipt_css = (REAL *)CHARM(malloc_aligned)(SIMD_MEMALIGN,
                                            size * sizeof(REAL));
    if (ipt_css == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }
    ipt_scs = (REAL *)CHARM(malloc_aligned)(SIMD_MEMALIGN,
                                            size * sizeof(REAL));
    if (ipt_scs == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }
    ipt_sss = (REAL *)CHARM(malloc_aligned)(SIMD_MEMALIGN,
                                            size * sizeof(REAL));
    if (ipt_sss == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }


    clt1v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                          sizeof(REAL));
    if (clt1v == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }
    clt2v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                          sizeof(REAL));
    if (clt2v == NULL)
    {
        FAILURE_priv = 1;
        goto FAILURE_1_parallel;
    }
    dcltv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                          sizeof(REAL));
    if (dcltv == NULL)
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
        if (CHARM(err_isempty)(err))
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);


        /* OK, and now all threads go to the "FAILURE_2_parallel" label to
         * deallocate all the memory that might be allocated before the
         * allocation failure. */
        goto FAILURE_2_parallel;
    }
    /* ..................................................................... */


    /* Lumped coefficients */
    REAL_SIMD lc00, lc01, lc10, lc11;


    /* Useful substitutions */
    unsigned int parity;
    unsigned int rem_m1_2, rem_m3_2;
    REAL dsigma;
    size_t row, ipv;
    size_t idx, idx2;
    REAL_SIMD tmp;


    /* Loop over the latitude cells */
#if CHARM_PARALLEL
#pragma omp for schedule(dynamic)
#endif
    for (size_t i = 0; i < SIMD_GET_MULTIPLE(cell_nlat); i += SIMD_SIZE)
    {
        /* Transformation of latitudes into co-latitudes */
        for (size_t v = 0; v < SIMD_SIZE; v++)
        {
            ipv = i + v;
            if (ipv < cell_nlat)
            {
                clt1v[v] = PI_2 - cell->latmax[ipv];
                clt2v[v] = PI_2 - cell->latmin[ipv];
            }
            else
            {
                clt1v[v] = clt2v[v] = PREC(0.0);
                continue;
            }


            dcltv[v] = clt2v[v] - clt1v[v];
        }


        /* Reset "fi" to zeros */
        memset(fi, 0, cell_nlon * SIMD_SIZE * sizeof(REAL));


        /* Reset "idx" to zero */
        idx = 0;


        /* Pre-computation of the trigonometric integrals */
        /* ----------------------------------------------------------------- */
        idx2 = 0;


        {
        REAL k1d, k3d;
        for (unsigned long k1 = 0; k1 <= nmax1; k1++)
        {
            k1d = (REAL)k1;
            for (unsigned long k3 = 0; k3 <= nmax3; k3++)
            {
                k3d = (REAL)k3;


                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    ipt_ccs[idx2] = CHARM(integ_ccs)(clt1v[v], dcltv[v], k1d,
                                                     k3d);
                    ipt_css[idx2] = CHARM(integ_css)(clt1v[v], dcltv[v], k1d,
                                                     k3d);
                    ipt_scs[idx2] = CHARM(integ_scs)(clt1v[v], dcltv[v], k1d,
                                                     k3d);
                    ipt_sss[idx2] = CHARM(integ_sss)(clt1v[v], dcltv[v], k1d,
                                                     k3d);


                    idx2++;
                }
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
            lc00 = lc01 = lc10 = lc11 = SET1_R(PREC(0.0));


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
                    tmp = LOAD_R(&ipt_ccs[idx2]);
                    lc00 = ADD_R(lc00, MUL_R(tmp, SET1_R(cnm1cnm3[idx])));
                    lc01 = ADD_R(lc01, MUL_R(tmp, SET1_R(cnm1snm3[idx])));
                    lc10 = ADD_R(lc10, MUL_R(tmp, SET1_R(snm1cnm3[idx])));
                    lc11 = ADD_R(lc11, MUL_R(tmp, SET1_R(snm1snm3[idx])));


                    idx++;
                    idx2 += SIMD_SIZE;
                }
                }
            }
            else if (parity == 1)
            /* (Even, Odd) */
            {
                for (unsigned long k1 = 0; k1 <= nmax1; k1++)
                {
                for (unsigned long k3 = 0; k3 <= nmax3; k3++)
                {
                    tmp = LOAD_R(&ipt_css[idx2]);
                    lc00 = ADD_R(lc00, MUL_R(tmp, SET1_R(cnm1cnm3[idx])));
                    lc01 = ADD_R(lc01, MUL_R(tmp, SET1_R(cnm1snm3[idx])));
                    lc10 = ADD_R(lc10, MUL_R(tmp, SET1_R(snm1cnm3[idx])));
                    lc11 = ADD_R(lc11, MUL_R(tmp, SET1_R(snm1snm3[idx])));


                    idx++;
                    idx2 += SIMD_SIZE;
                }
                }
            }
            else if (parity == 2)
            /* (Odd, Even) */
            {
                for (unsigned long k1 = 0; k1 <= nmax1; k1++)
                {
                for (unsigned long k3 = 0; k3 <= nmax3; k3++)
                {
                    tmp = LOAD_R(&ipt_scs[idx2]);
                    lc00 = ADD_R(lc00, MUL_R(tmp, SET1_R(cnm1cnm3[idx])));
                    lc01 = ADD_R(lc01, MUL_R(tmp, SET1_R(cnm1snm3[idx])));
                    lc10 = ADD_R(lc10, MUL_R(tmp, SET1_R(snm1cnm3[idx])));
                    lc11 = ADD_R(lc11, MUL_R(tmp, SET1_R(snm1snm3[idx])));


                    idx++;
                    idx2 += SIMD_SIZE;
                }
                }
            }
            else /*if (parity == 3) */
            /* (Odd, Odd) */
            {
                for (unsigned long k1 = 0; k1 <= nmax1; k1++)
                {
                for (unsigned long k3 = 0; k3 <= nmax3; k3++)
                {
                    tmp = LOAD_R(&ipt_sss[idx2]);
                    lc00 = ADD_R(lc00, MUL_R(tmp, SET1_R(cnm1cnm3[idx])));
                    lc01 = ADD_R(lc01, MUL_R(tmp, SET1_R(cnm1snm3[idx])));
                    lc10 = ADD_R(lc10, MUL_R(tmp, SET1_R(snm1cnm3[idx])));
                    lc11 = ADD_R(lc11, MUL_R(tmp, SET1_R(snm1snm3[idx])));


                    idx++;
                    idx2 += SIMD_SIZE;
                }
                }
            }
            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


            /* Synthesis for the "i"th latitude parallel */
            CHARM(shs_cell_isurf_lr)(lon0, dlon, cell_nlon,
                                     lc00, lc01, lc10, lc11, m1, m3,
                                     fi);
        }  /* End of the loop over "m3" */
        }  /* End of the loop over "m1" */


        /* Final synthesis */
        /* --------------------------------------------------------- */
        for (size_t v = 0; v < SIMD_SIZE; v++)
        {
            ipv = i + v;
            if (ipv >= cell_nlat)
                continue;


            row = ipv * cell_nlon;
            dcltv[v] = COS(clt1v[v]) - COS(clt2v[v]);


            for (size_t j = 0; j < cell_nlon; j++)
            {
                /* Compute the area of the cells on the unit sphere */
                dsigma = dcltv[v] * DELTAlon[j];


                /* Final synthesis */
                f[row + j] = (mur / dsigma) * fi[j * SIMD_SIZE + v];
            }
        }
        /* --------------------------------------------------------- */


    }  /* End of the loop over "i" */


    /* Free the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE_2_parallel:
    CHARM(free_aligned)(fi);
    CHARM(free_aligned)(ipt_ccs);
    CHARM(free_aligned)(ipt_css);
    CHARM(free_aligned)(ipt_scs);
    CHARM(free_aligned)(ipt_sss);
    CHARM(free_aligned)(clt1v);
    CHARM(free_aligned)(clt2v);
    CHARM(free_aligned)(dcltv);
    /* --------------------------------------------------------------------- */


    }
    /* --------------------------------------------------------------------- */






    /* Free the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if ((FAILURE_glob != 0) && CHARM(err_isempty)(err))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(DELTAlon);
    CHARM(free_aligned)(cnm1cnm3);
    CHARM(free_aligned)(cnm1snm3);
    CHARM(free_aligned)(snm1cnm3);
    CHARM(free_aligned)(snm1snm3);
    /* --------------------------------------------------------------------- */





    return;
}
