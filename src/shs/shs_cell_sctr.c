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
#include "shs_cell_kernel.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_gm_hm.h"
#include "../leg/leg_pol_en_fn.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../err/err_set.h"
#include "shs_rpows.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_sctr)(const CHARM(crd) *cell, const CHARM(shc) *shcs,
                          unsigned long nmax, REAL *f, CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    int FAILURE_glob = 0;
    REAL *r = NULL;
    REAL *ri = NULL;
    REAL *dm = NULL;
    REAL *en = NULL;
    REAL *fn = NULL;
    REAL *gm = NULL;
    REAL *hm = NULL;
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


    /* Coefficients "dm" for recurrence relations to compute fully-normalized
     * Legendre functions */
    dm = (REAL *)calloc(nmax + 1, sizeof(REAL));
    if (dm == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    CHARM(leg_func_dm)(nmax, r, ri, dm);


    /* Coefficients "en" and "fn" for recurrence relations to compute
     * *un-normalized* Legendre polynomials. (Note that to compute integrals of
     * fully-normalized Legendre polynomials, we need un-normalized Legendre
     * polynomials; see, e.g., Eq. B.24 of Jekeli et al. 2007) */
    en = (REAL *)calloc(nmax + 2, sizeof(REAL));
    if (en == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    fn = (REAL *)calloc(nmax + 2, sizeof(REAL));
    if (fn == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    CHARM(leg_pol_en_fn)(nmax + 1, en, fn);
    /* --------------------------------------------------------------------- */






    /* Initializations for recurrence relations to compute integrals of
     * fully-normalized Legendre functions */
    /* --------------------------------------------------------------------- */
    gm = (REAL *)calloc(nmax + 1, sizeof(REAL));
    if (gm == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    hm = (REAL *)calloc(nmax + 1, sizeof(REAL));
    if (hm == NULL)
    {
        FAILURE_glob = 1;
        goto FAILURE;
    }
    CHARM(leg_func_gm_hm)(nmax, r, ri, gm, hm);
    /* --------------------------------------------------------------------- */






    /* Get the number of evaluation cells */
    /* --------------------------------------------------------------------- */
    /* We have already checked that "cell->nlat == cell->nlon", so any of the
     * two can be used to get the number of evaluation points. */
    size_t ncells = cell->nlat;
    /* --------------------------------------------------------------------- */






    /* Loop over grid latitudes */
    /* --------------------------------------------------------------------- */
    REAL mur = shcs->mu / shcs->r;


#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(f, shcs, nmax, cell, dm, en, fn, gm, hm, r, ri) \
shared(ncells, FAILURE_glob, mur, err)
#endif
    {
        /* ................................................................. */
        /* An indicator for failed memory initializations on each thread,
         * a private variable. */
        int FAILURE_priv = 0;


        int   *ips1   = NULL;
        int   *ips2   = NULL;
        REAL *ps1   = NULL;
        REAL *ps2   = NULL;
        REAL *anm   = NULL;
        REAL *bnm   = NULL;
        REAL *rpows = NULL;


        ips1 =    (int *)calloc(nmax, sizeof(int));
        if (ips1 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        ips2 =    (int *)calloc(nmax, sizeof(int));
        if (ips2 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        ps1  = (REAL *)calloc(nmax, sizeof(REAL));
        if (ps1 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        ps2  = (REAL *)calloc(nmax, sizeof(REAL));
        if (ps2 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        anm  = (REAL *)calloc(nmax + 1, sizeof(REAL));
        if (anm == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        bnm  = (REAL *)calloc(nmax + 1, sizeof(REAL));
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
            /* At this point, "FAILURE_priv" is equal to 1 if any
             * allocation failed on a particular thread.  Now, we can add
             * "FAILURE_priv" from all the threads to get the total number
             * threads on which an allocation failure occurred (the
             * "FAILURE_glob" variable).  Note that this code block is
             * executed with one thread at a time only, so "FAILURE_glob"
             * can safely be overwritten.   */
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


        REAL latmini, latmaxi;
        REAL lon1i, lon2i, dlon;
        REAL t1, u1;
        REAL t2, u2;

        REAL a, b, a2, b2;
        a = b = a2 = b2 = ADDP(0.0);
        REAL imm0, imm1, imm2;
        REAL lontmp, m2, m2sm2dl;
        REAL dsigma;
        REAL fi = ADDP(0.0);


#if CHARM_PARALLEL
#pragma omp for
#endif
        for (size_t i = 0; i < ncells; i++)
        {
            /* Pre-compute the powers of "shcs->r / pnt->r[i]" */
            /* ------------------------------------------------------------- */
            CHARM(shs_rpows)(shcs->r, cell->r[i], rpows, nmax);
            /* ------------------------------------------------------------- */


            /* Prepare arrays for sectorial Legendre functions */
            /* ------------------------------------------------------------- */
            latmini = cell->lat[2 * i];
            latmaxi = cell->lat[2 * i + 1];


            t1    = SIN(latmini);
            u1    = COS(latmini);
            CHARM(leg_func_prepare)(u1, ps1, ips1, dm, nmax);


            t2    = SIN(latmaxi);
            u2    = COS(latmaxi);
            CHARM(leg_func_prepare)(u2, ps2, ips2, dm, nmax);
            /* ------------------------------------------------------------- */


            /* The "fi" variable represents the synthesized quantity "f" for
             * the "i"th latitude and longitude. Therefore, it needs to be
             * reinitialized to zero. */
            fi  = ADDP(0.0);


            /* Some useful substitutions */
            lon1i = cell->lon[2 * i];
            lon2i = cell->lon[2 * i + 1];
            dlon  = lon2i - lon1i;


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {


                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Summation over harmonic degrees */
                CHARM(shs_cell_kernel)(nmax, m, shcs, anm, bnm,
                                       latmini, latmaxi,
                                       t1, t2,
                                       u1, u2,
                                       ps1, ps2,
                                       ips1, ips2,
                                       &imm0, &imm1, &imm2,
                                       en, fn, gm, hm, ri, rpows, NULL,
                                       0,
                                       &a, &b, &a2, &b2);


                /* The longitudinal part of the synthesis */
                /* --------------------------------------------------------- */
                if (m == 0)
                    m2sm2dl = dlon;
                else
                {
                    m2 = ADDP(2.0) / (REAL)m;
                    m2sm2dl = m2 * SIN(dlon / m2);
                }

                a *= m2sm2dl;
                b *= m2sm2dl; /* Remember that "b = 0.0" for "m == 0.0", so it
                               * can safely be multiplied by "m2sm2dl = dlon"
                               * */


                lontmp = (REAL)m * (lon1i + lon2i) / ADDP(2.0);
                fi    += a * COS(lontmp) + b * SIN(lontmp);
                /* --------------------------------------------------------- */


            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            /* Synthesis */
            /* ------------------------------------------------------------- */
            /* Area of the cell on the unit sphere */
            dsigma = (SIN(latmaxi) - SIN(latmini)) * dlon;


            f[i] = mur * fi / dsigma;
            /* ------------------------------------------------------------- */


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


FAILURE_2_parallel:
        free(ips1); free(ps1);
        free(ips2); free(ps2);
        free(anm); free(bnm);
        free(rpows);


    } /* End of "#pragma omp parallel" */
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if (FAILURE_glob != 0)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(r); free(ri);
    free(dm); free(en); free(fn); free(gm); free(hm);
    /* --------------------------------------------------------------------- */






    return;
}
