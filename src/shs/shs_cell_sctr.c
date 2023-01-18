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
#include "shs_rpows.h"
#include "shs_sctr_mulc.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_gm_hm.h"
#include "../leg/leg_pol_en_fn.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../crd/crd_check_cells.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_sctr)(const CHARM(cell) *cell, const CHARM(shc) *shcs,
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






    /* Check cell boundaries */
    /* ..................................................................... */
    CHARM(crd_check_cells)(cell, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* ..................................................................... */






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
    size_t ncell = cell->nlat;
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    REAL_SIMD mur = SET1_R(shcs->mu / shcs->r);


    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(f, shcs, nmax, cell, dm, en, fn, gm, hm, r, ri, pt) \
shared(ncell, FAILURE_glob, mur, err)
#endif
    {
        /* ................................................................. */
        /* An indicator for failed memory initializations on each thread,
         * a private variable. */
        int FAILURE_priv = 0;


        int  *ips1    = NULL;
        int  *ips2    = NULL;
        REAL *ps1     = NULL;
        REAL *ps2     = NULL;
        REAL *t1v     = NULL;
        REAL *u1v     = NULL;
        REAL *t2v     = NULL;
        REAL *u2v     = NULL;
        REAL *latminv = NULL;
        REAL *latmaxv = NULL;
        REAL *lonminv = NULL;
        REAL *lonmaxv = NULL;
        REAL *clonv   = NULL;
        REAL *slonv   = NULL;
        REAL *dlonv   = NULL;
        REAL *tmpv    = NULL;
        REAL *anm     = NULL;
        REAL *bnm     = NULL;
        REAL *rpows   = NULL;


        ips1 = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nmax * SIMD_SIZE,
                                            sizeof(int));
        if (ips1 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        ips2 = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nmax * SIMD_SIZE,
                                            sizeof(int));
        if (ips2 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        ps1 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nmax * SIMD_SIZE,
                                             sizeof(REAL));
        if (ps1 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        ps2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nmax * SIMD_SIZE,
                                             sizeof(REAL));
        if (ps2 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        t1v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (t1v == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        u1v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (u1v == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        t2v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (t2v == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        u2v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (u2v == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        latminv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        if (latminv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        latmaxv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        if (latmaxv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        lonminv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        if (lonminv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        lonmaxv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        if (lonmaxv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        clonv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                              sizeof(REAL));
        if (clonv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        slonv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                              sizeof(REAL));
        if (slonv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        dlonv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                              sizeof(REAL));
        if (dlonv == NULL)
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
            /* At this point, "FAILURE_priv" is equal to 1 if any
             * allocation failed on a particular thread.  Now, we can add
             * "FAILURE_priv" from all the threads to get the total number
             * threads on which an allocation failure occurred (the
             * "FAILURE_glob" variable).  Note that this code block is
             * executed with one thread at a time only, so "FAILURE_glob" can
             * safely be overwritten. */
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


        REAL_SIMD latmin, latmax, dlon, t1, u1, t2, u2;

        REAL_SIMD a, b, a2, b2;
        a = b = a2 = b2 = SET_ZERO_R;
        REAL_SIMD imm0, imm1, imm2;
        REAL_SIMD m2sm2dl, fi, clon, slon;
        REAL m2, m2i, lontmp;
        REAL_SIMD dsigma;
        REAL_SIMD tmp = SET_ZERO_R;
;
        size_t ipv;


#if CHARM_PARALLEL
#pragma omp for schedule(dynamic)
#endif
        for (size_t i = 0; i < SIMD_GET_MULTIPLE(ncell); i += SIMD_SIZE)
        {
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                ipv = i + v;
                if (ipv < ncell)
                {
                    latmaxv[v] = cell->latmax[ipv];
                    latminv[v] = cell->latmin[ipv];


                    t1v[v] = SIN(latminv[v]);
                    u1v[v] = COS(latminv[v]);


                    t2v[v] = SIN(latmaxv[v]);
                    u2v[v] = COS(latmaxv[v]);


                    lonminv[v] = cell->lonmin[ipv];
                    lonmaxv[v] = cell->lonmax[ipv];
                    dlonv[v]   = lonmaxv[v] - lonminv[v];
                }
                else
                {
                    latminv[v] = latmaxv[v] = t1v[v] = u1v[v] = 
                        t2v[v] = u2v[v] = lonminv[v] = lonmaxv[v] = 
                        dlonv[v] = PREC(0.0);
                    continue;
                }


                /* Pre-compute the powers of "shcs->r / pnt->r[ipv]" */
                CHARM(shs_rpows)(v, shcs->r, cell->r[ipv], rpows, nmax);
            }


            t1     = LOAD_R(&t1v[0]);
            t2     = LOAD_R(&t2v[0]);
            u1     = LOAD_R(&u1v[0]);
            u2     = LOAD_R(&u2v[0]);
            latmin = LOAD_R(&latminv[0]);
            latmax = LOAD_R(&latmaxv[0]);
            dlon   = LOAD_R(&dlonv[0]);
            fi     = SET_ZERO_R;


            CHARM(leg_func_prepare)(u1v, ps1, ips1, dm, nmax);
            CHARM(leg_func_prepare)(u2v, ps2, ips2, dm, nmax);


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Apply polar optimization if asked to do so.  Since "u1
                 * = sin(latmin)" and "u2 = sin(latmax)", it is sufficient to
                 * check "u1" only.  In other words, if the polar optimization
                 * can be applied for "u1", it can surely be applied for
                 * "u2". */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, u1, pt))
                    continue;


                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Summation over harmonic degrees */
                CHARM(shs_cell_kernel)(nmax, m, shcs, anm, bnm,
                                       latmin, latmax,
                                       t1, t2,
                                       u1, u2,
                                       ps1, ps2,
                                       ips1, ips2,
                                       &imm0, &imm1, &imm2,
                                       en, fn, gm, hm, ri, rpows, NULL,
                                       SET_ZERO_R,
                                       &a, &b, &a2, &b2);


                /* The longitudinal part of the synthesis */
                /* --------------------------------------------------------- */
                if (m == 0)
                    m2sm2dl = dlon;
                else
                {
                    m2 = PREC(2.0) / (REAL)m;
                    for (size_t v = 0; v < SIMD_SIZE; v++)
                        tmpv[v] = m2 * SIN(dlonv[v] / m2);
                    m2sm2dl = LOAD_R(&tmpv[0]);
                }


                a = MUL_R(a, m2sm2dl);
                b = MUL_R(b, m2sm2dl); /* Remember that "b = 0.0" for "m ==
                                        * 0.0", so it can safely be multiplied
                                        * by "m2sm2dl = dlon" */


                m2i = (REAL)m / PREC(2.0);
                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    lontmp = m2i * (lonminv[v] + lonmaxv[v]);
                    clonv[v] = COS(lontmp);
                    slonv[v] = SIN(lontmp);
                }
                clon = LOAD_R(&clonv[0]);
                slon = LOAD_R(&slonv[0]);
                fi = ADD_R(fi, ADD_R(MUL_R(a, clon), MUL_R(b, slon)));
                /* --------------------------------------------------------- */


            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            /* Final part of the synthesis */
            /* ------------------------------------------------------------- */
            /* Area of the cell on the unit sphere */
            dsigma = MUL_R(SUB_R(t2, t1), dlon);


            CHARM(shs_sctr_mulc)(i, ncell, mur, tmp, tmpv, DIV_R(fi, dsigma),
                                 f);
            /* ------------------------------------------------------------- */


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


FAILURE_2_parallel:
        CHARM(free_aligned)(ips1);    CHARM(free_aligned)(ps1);
        CHARM(free_aligned)(ips2);    CHARM(free_aligned)(ps2);
        CHARM(free_aligned)(t1v);     CHARM(free_aligned)(u1v);
        CHARM(free_aligned)(t2v);     CHARM(free_aligned)(u2v);
        CHARM(free_aligned)(latminv); CHARM(free_aligned)(latmaxv);
        CHARM(free_aligned)(lonminv); CHARM(free_aligned)(lonmaxv);
        CHARM(free_aligned)(dlonv);   CHARM(free_aligned)(tmpv);
        CHARM(free_aligned)(clonv);   CHARM(free_aligned)(slonv);
        CHARM(free_aligned)(rpows);
        free(anm); free(bnm);


    } /* End of "#pragma omp parallel" */
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if ((FAILURE_glob != 0) && CHARM(err_isempty)(err))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(r); free(ri);
    free(dm); free(en); free(fn); free(gm); free(hm);
    /* --------------------------------------------------------------------- */






    return;
}
