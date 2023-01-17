/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
#include "shs_point_kernel.h"
#include "shs_rpows.h"
#include "shs_sctr_mulc.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#include "../err/err_set.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_point_sctr)(const CHARM(point) *pnt, const CHARM(shc) *shcs,
                           unsigned long nmax, REAL *f, CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    int FAILURE_glob = 0;
    REAL *r  = NULL;
    REAL *ri = NULL;
    REAL *dm = NULL;
    /* --------------------------------------------------------------------- */






    /* Initializations for recurrence relations to compute Legendre
     * functions. */
    /* --------------------------------------------------------------------- */
    /* Prepare some variables to compute coefficients the "anm" and "bnm"
     * coefficients for the Legendre functions recurrences */
    r = (REAL *)calloc(2 * nmax + 4, sizeof(REAL));
    if (r == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    ri = (REAL *)calloc(2 * nmax + 4, sizeof(REAL));
    if (ri == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    CHARM(leg_func_r_ri)(nmax, r, ri);


    /* "dm" coefficients for sectorial Legendre functions */
    dm = (REAL *)calloc(nmax + 1, sizeof(REAL));
    if (dm == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    CHARM(leg_func_dm)(nmax, r, ri, dm);
    /* --------------------------------------------------------------------- */






    /* Get the number of evaluation points */
    /* --------------------------------------------------------------------- */
    /* We have already checked that "pnt->nlat == pnt->nlon", so any of the two
     * can be used to get the number of evaluation points. */
    size_t npnt = pnt->nlat;
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    REAL_SIMD mur = SET1_R(shcs->mu / shcs->r);


    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(f, shcs, nmax, pnt, npnt, dm, r, ri, FAILURE_glob, mur, err, pt)
#endif
    {
        /* ................................................................. */
        /* An indicator for failed memory initializations on each thread,
         * a private variable. */
        int FAILURE_priv = 0;


        int  *ips     = NULL;
        REAL *ps      = NULL;
        REAL *tv      = NULL;
        REAL *uv      = NULL;
        REAL *lonv    = NULL;
        REAL *clonimv = NULL;
        REAL *slonimv = NULL;
        REAL *tmpv    = NULL;
        REAL *anm     = NULL;
        REAL *bnm     = NULL;
        REAL *rpows   = NULL;


        ips = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nmax * SIMD_SIZE,
                                           sizeof(int));
        if (ips == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        ps = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nmax * SIMD_SIZE,
                                           sizeof(REAL));
        if (ps == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        tv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                           sizeof(REAL));
        if (tv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        uv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                           sizeof(REAL));
        if (uv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        lonv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                             sizeof(REAL));
        if (lonv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        clonimv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        if (clonimv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        slonimv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        if (slonimv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        tmpv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                             sizeof(REAL));
        if (tmpv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        anm = (REAL *)calloc(nmax + 1, sizeof(REAL));
        if (anm == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        bnm = (REAL *)calloc(nmax + 1, sizeof(REAL));
        if (bnm == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        rpows = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                              (nmax + 1) * SIMD_SIZE,
                                              sizeof(REAL));
        if (rpows == NULL)
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
             * that this code block is executed with one thread at a time only,
             * so "FAILURE_glob" can safely be overwritten.   */
            FAILURE_glob += FAILURE_priv;
        }


        /* Now we have to wait until all the threads get here. */
#if CHARM_PARALLEL
#pragma omp barrier
#endif
        /* OK, now let's check on each thread whether there is at least one
         * failed memory allocation among the threads. */
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


            /* OK, and now all threads go to the "FAILURE_2_parallel" label
             * to deallocate all the memory that might be allocated before the
             * allocation failure. */
            goto FAILURE_2_parallel;
        }
        /* ................................................................. */


        REAL_SIMD t, u, fi;
        REAL_SIMD a, b, a2, b2;
        a = b = a2 = b2 = SET_ZERO_R;
        REAL_SIMD clonim, slonim;
        REAL_SIMD tmp = SET_ZERO_R;
        REAL lontmp;
        size_t ipv;


#if CHARM_PARALLEL
#pragma omp for
#endif
        for (size_t i = 0; i < SIMD_GET_MULTIPLE(npnt); i += SIMD_SIZE)
        {
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                ipv = i + v;
                if (ipv < npnt)
                {
                    tv[v]   = SIN(pnt->lat[ipv]);
                    uv[v]   = COS(pnt->lat[ipv]);
                    lonv[v] = pnt->lon[ipv];
                }
                else
                {
                    tv[v] = uv[v] = lonv[v] = PREC(0.0);
                    continue;
                }


                /* Pre-compute the powers of "shcs->r / pnt->r[ipv]" */
                CHARM(shs_rpows)(v, shcs->r, pnt->r[ipv], rpows, nmax);
            }


            t  = LOAD_R(&tv[0]);
            u  = LOAD_R(&uv[0]);
            fi = SET_ZERO_R;


            /* Prepare arrays for sectorial Legendre functions */
            CHARM(leg_func_prepare)(uv, ps, ips, dm, nmax);


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Apply polar optimization if asked to do so */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, u, pt))
                    continue;


                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Summation over harmonic degrees.  Note that the symmetry of
                 * Legendre functions cannot be utilized with scattered points,
                 * hence "SET_ZERO_R".  The output variables "a2" and "b2" are
                 * not used with scattered points. */
                CHARM(shs_point_kernel)(nmax, m, shcs, anm, bnm,
                                        t, ps, ips, rpows, NULL,
                                        SET_ZERO_R, &a, &b, &a2, &b2);


                /* The longitudinal part of the synthesis */
                /* ......................................................... */
                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    lontmp = (REAL)m * lonv[v];
                    clonimv[v] = COS(lontmp);
                    slonimv[v] = SIN(lontmp);
                }
                clonim = LOAD_R(&clonimv[0]);
                slonim = LOAD_R(&slonimv[0]);


                fi = ADD_R(fi, ADD_R(MUL_R(a, clonim), MUL_R(b, slonim)));
                /* ......................................................... */


            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            /* Final part of the synthesis */
            CHARM(shs_sctr_mulc)(i, npnt, mur, tmp, tmpv, fi, f);


        } /* End of the loop over the evaluation points */
        /* ----------------------------------------------------------------- */


        /* Free the heap memory */
        /* ----------------------------------------------------------------- */
FAILURE_2_parallel:
        CHARM(free_aligned)(ips);      CHARM(free_aligned)(ps);
        CHARM(free_aligned)(tv);       CHARM(free_aligned)(uv);
        CHARM(free_aligned)(rpows);    CHARM(free_aligned)(lonv);
        CHARM(free_aligned)(clonimv);  CHARM(free_aligned)(slonimv);
        CHARM(free_aligned)(tmpv);
        free(anm); free(bnm);
        /* ----------------------------------------------------------------- */


    } /* End of "#pragma omp parallel" */
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if ((FAILURE_glob != 0) && CHARM(err_isempty)(err))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(r); free(ri); free(dm);
    /* --------------------------------------------------------------------- */






    return;
}
