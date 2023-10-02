/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "shs_point_kernel.h"
#include "shs_sctr_mulc.h"
#include "shs_r_eq_rref.h"
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






    /* Check whether all values of "pnt->r" are equal to "shcs->r".  If true,
     * a faster code can be used inside "shs_point_kernel".  */
    /* --------------------------------------------------------------------- */
    _Bool r_eq_rref = CHARM(shs_r_eq_rref)(pnt, shcs);
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    REAL_SIMD mur = SET1_R(shcs->mu / shcs->r);


    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


    /* Radius of the reference sphere that is associated with the spherical
     * harmonic coefficients */
    REAL_SIMD rref = SET1_R(shcs->r);


#if CHARM_OPENMP
#pragma omp parallel default(none) \
shared(f, shcs, nmax, pnt, npnt, dm, r, ri, FAILURE_glob, mur, err, pt, rref) \
shared(r_eq_rref)
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
        REAL *pnt_rv  = NULL;
        REAL *tmpv    = NULL;
        REAL *anm     = NULL;
        REAL *bnm     = NULL;


        ips = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           nmax * SIMD_SIZE * SIMD_BLOCK,
                                           sizeof(int));
        if (ips == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        ps = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           nmax * SIMD_SIZE * SIMD_BLOCK,
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
        lonv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                             SIMD_SIZE * SIMD_BLOCK,
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
        pnt_rv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                               sizeof(REAL));
        if (pnt_rv == NULL)
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


FAILURE_1_parallel:
#if CHARM_OPENMP
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
#if CHARM_OPENMP
#pragma omp barrier
#endif
        /* OK, now let's check on each thread whether there is at least one
         * failed memory allocation among the threads. */
        if (FAILURE_glob > 0)
        {
            /* Ooops, there was indeed a memory allocation failure.  So let the
             * master thread write down the error to the "err" variable. */
#if CHARM_OPENMP
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


        size_t ipv, l;


        REAL_SIMD t[SIMD_BLOCK], u[SIMD_BLOCK], fi[SIMD_BLOCK];
        REAL_SIMD pnt_r[SIMD_BLOCK];
        REAL_SIMD ratio[SIMD_BLOCK], ratiom[SIMD_BLOCK];
        REAL_SIMD a[SIMD_BLOCK], b[SIMD_BLOCK], a2[SIMD_BLOCK], b2[SIMD_BLOCK];
        for (l = 0; l < SIMD_BLOCK; l++)
            a[l] = b[l] = a2[l] = b2[l] = ratio[l] = ratiom[l] = SET_ZERO_R;
        REAL_SIMD zeros;
        REAL_SIMD clonim, slonim;
        REAL_SIMD tmp = SET_ZERO_R;
        REAL lontmp, m_real;


#if CHARM_OPENMP
#pragma omp for schedule(dynamic)
#endif
        for (size_t i = 0; i < SIMD_MULTIPLE(npnt, SIMD_SIZE * SIMD_BLOCK);
             i += SIMD_SIZE * SIMD_BLOCK)
        {
            for (l = 0; l < SIMD_BLOCK; l++)
            {
                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    ipv = i + l * SIMD_SIZE + v;
                    if (ipv < npnt)
                    {
                        tv[v]     = SIN(pnt->lat[ipv]);
                        uv[v]     = COS(pnt->lat[ipv]);
                        lonv[l * SIMD_SIZE + v]   = pnt->lon[ipv];
                        pnt_rv[v] = pnt->r[ipv];
                    }
                    else
                    {
                        tv[v] = uv[v] = lonv[l * SIMD_SIZE + v] = pnt_rv[v] =
                            PREC(0.0);
                        continue;
                    }
                }


                t[l]     = LOAD_R(&tv[0]);
                u[l]     = LOAD_R(&uv[0]);
                pnt_r[l] = LOAD_R(&pnt_rv[0]);
                fi[l]    = SET_ZERO_R;


                ratio[l]  = DIV_R(rref, pnt_r[l]);
                ratiom[l] = ratio[l];


                /* Prepare arrays for sectorial Legendre functions */
                CHARM(leg_func_prepare)(uv, ps + l * SIMD_SIZE * nmax,
                                        ips + l * SIMD_SIZE * nmax, dm, nmax);
            }


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Apply polar optimization if asked to do so */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, &u[0],
                                                         SIMD_BLOCK, pt))
                    goto UPDATE_RATIOS;


                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Summation over harmonic degrees.  Note that the symmetry of
                 * Legendre functions cannot be utilized with scattered points,
                 * hence "SET_ZERO_R".  The output variables "a2" and "b2" are
                 * not used with scattered points. */
                CHARM(shs_point_kernel)(nmax, m, shcs, r_eq_rref, anm, bnm,
                                        &t[0], ps, ips,
                                        &ratio[0], &zeros,
                                        &ratiom[0], &zeros,
                                        &zeros,
                                        &a[0], &b[0], &a2[0], &b2[0]);


                /* The longitudinal part of the synthesis */
                /* ......................................................... */
                m_real = (REAL)m;


                for (l = 0; l < SIMD_BLOCK; l++)
                {
                    for (size_t v = 0; v < SIMD_SIZE; v++)
                    {
                        lontmp = m_real * lonv[l * SIMD_SIZE + v];
                        clonimv[v] = COS(lontmp);
                        slonimv[v] = SIN(lontmp);
                    }
                    clonim = LOAD_R(&clonimv[0]);
                    slonim = LOAD_R(&slonimv[0]);


                    fi[l] = ADD_R(fi[l],
                                  ADD_R(MUL_R(a[l], clonim),
                                        MUL_R(b[l], slonim)));
                }
                /* ......................................................... */


UPDATE_RATIOS:
                for (l = 0; l < SIMD_BLOCK; l++)
                    ratiom[l] = MUL_R(ratiom[l], ratio[l]);


            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            /* Final part of the synthesis */
            CHARM(shs_sctr_mulc)(i, npnt, pnt->type, mur, tmp, tmpv, &fi[0],
                                 f);


        } /* End of the loop over the evaluation points */
        /* ----------------------------------------------------------------- */


        /* Free the heap memory */
        /* ----------------------------------------------------------------- */
FAILURE_2_parallel:
        CHARM(free_aligned)(ips);      CHARM(free_aligned)(ps);
        CHARM(free_aligned)(tv);       CHARM(free_aligned)(uv);
        CHARM(free_aligned)(pnt_rv);   CHARM(free_aligned)(lonv);
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
