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
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../err/err_set.h"
#include "shs_rpows.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_point_sctr)(const CHARM(crd) *pnt, const CHARM(shc) *shcs,
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
    REAL mur = shcs->mu / shcs->r;


#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(f, shcs, nmax, pnt, npnt, dm, r, ri, FAILURE_glob, mur, err)
#endif
    {
        /* ................................................................. */
        /* An indicator for failed memory initializations on each thread,
         * a private variable. */
        int FAILURE_priv = 0;


        int    *ips   = NULL;
        REAL  *ps   = NULL;
        REAL *anm   = NULL;
        REAL *bnm   = NULL;
        REAL *rpows = NULL;


        ips = (int *)calloc(nmax, sizeof(int));
        if (ips == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        ps  = (REAL *)calloc(nmax, sizeof(REAL));
        if (ps == NULL)
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


        rpows = (REAL *)calloc(nmax + 1, sizeof(REAL));
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
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                               CHARM_ERR_MALLOC_FAILURE);
            }


            /* OK, and now all threads go to the "FAILURE_2_parallel" label
             * to deallocate all the memory that might be allocated before the
             * allocation failure. */
            goto FAILURE_2_parallel;
        }
        /* ................................................................. */


        REAL fi;
        REAL t, u;

        REAL a, b, a2, b2;
        a = b = a2 = b2 = ADDP(0.0);
        REAL lonim;


#if CHARM_PARALLEL
#pragma omp for
#endif
        for (size_t i = 0; i < npnt; i++)
        {
            /* Prepare arrays for sectorial Legendre functions */
            /* ------------------------------------------------------------- */
            t  = SIN(pnt->lat[i]);
            u  = COS(pnt->lat[i]);
            CHARM(leg_func_prepare)(u, ps, ips, dm, nmax);
            /* ------------------------------------------------------------- */


            /* Pre-compute the powers of "shcs->r / pnt->r[i]" */
            /* ------------------------------------------------------------- */
            CHARM(shs_rpows)(shcs->r, pnt->r[i], rpows, nmax);
            /* ------------------------------------------------------------- */


            /* The "fi" variable represents the synthesized quantity "f" for
             * the "i"th latitude and longitude. Therefore, it needs to be
             * reinitialized to zero. */
            fi = ADDP(0.0);


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Summation over harmonic degrees.  Note that the symmetry of
                 * Legendre functions cannot be utilized with scattered points,
                 * hence "symmi = 0".  The output variables "a2" and "b2" are
                 * not used with scattered points. */
                CHARM(shs_point_kernel)(nmax, m, shcs, anm, bnm,
                                        t, ps, ips, rpows, NULL,
                                        0,
                                        &a, &b, &a2, &b2);


                /* The longitudinal part of the synthesis */
                lonim = (REAL)m * pnt->lon[i];
                fi  += a * COS(lonim) + b * SIN(lonim);


            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            /* Save the synthesized "i"th value of "f" */
            f[i] = mur * fi;


        } /* End of the loop over the evaluation points */
        /* ----------------------------------------------------------------- */


        /* Free the heap memory */
        /* ----------------------------------------------------------------- */
FAILURE_2_parallel:
        free(ips); free(ps);
        free(anm); free(bnm);
        free(rpows);
        /* ----------------------------------------------------------------- */


    } /* End of "#pragma omp parallel" */
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if (FAILURE_glob != 0)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(r); free(ri); free(dm);
    /* --------------------------------------------------------------------- */






    return;
}
