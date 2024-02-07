/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_enm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
#include "shs_point_kernels.h"
#include "shs_sctr_mulc.h"
#include "shs_r_eq_rref.h"
#include "shs_lc_struct.h"
#include "shs_lc_init.h"
#include "shs_get_mur_dorder_npar.h"
#include "shs_point_gradn.h"
#include "shs_max_npar.h"
#include "shs_point_sctr.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
#undef LC_CS
#define LC_CS(x)                                                              \
        fi[idx] = ADD_R(fi[idx], ADD_R(MUL_R(lc.CAT(a, x)[l], clonim),        \
                                       MUL_R(lc.CAT(b, x)[l], slonim)));
/* ------------------------------------------------------------------------- */






void CHARM(shs_point_sctr)(const CHARM(point) *pnt,
                           const CHARM(shc) *shcs,
                           unsigned long nmax,
                           int dr,
                           int dlat,
                           int dlon,
                           REAL **f,
                           CHARM(err) *err)
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






    /* Get some constants */
    /* --------------------------------------------------------------------- */
    REAL mur;  /* "(shcs->mu / shcs->r)^dorder" */
    unsigned dorder;  /* "0" for potential,
                         "1" for first-order derivatives,
                         "2" for second-order derivatives */
    size_t npar; /* Number quantities to be synthesized:
                    "3" for first-order derivatives,
                    "6" for second-order derivatives,
                    "1" otherwise. */
    CHARM(shs_get_mur_dorder_npar)(shcs, dr, dlat, dlon, &mur, &dorder, &npar,
                                   err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    int grad;
    if ((dr == GRAD_1) && (dlat == GRAD_1) && (dlon == GRAD_1))
        grad = 1;
    else if ((dr == GRAD_2) && (dlat == GRAD_2) && (dlon == GRAD_2))
        grad = 2;
    else
        grad = 0;


    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


    /* Radius of the reference sphere that is associated with the spherical
     * harmonic coefficients */
    REAL_SIMD rref = SET1_R(shcs->r);


#if CHARM_OPENMP
#pragma omp parallel default(none) \
shared(f, shcs, nmax, pnt, npnt, dm, r, ri, FAILURE_glob, mur, err, pt, rref) \
shared(r_eq_rref, dr, dlat, dlon, npar, grad, dorder)
#endif
    {
        REAL_SIMD fi[SHS_MAX_NPAR * SIMD_BLOCK];


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
        REAL *enm     = NULL;


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
        if (dorder > 0)
        {
            enm = (REAL *)calloc(nmax + 1, sizeof(REAL));
            if (enm == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
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


        size_t ipv, l, idx;


        REAL_SIMD t[SIMD_BLOCK], u[SIMD_BLOCK];
        REAL_SIMD pnt_r[SIMD_BLOCK];
        REAL_SIMD ratio[SIMD_BLOCK], ratiom[SIMD_BLOCK];
        for (l = 0; l < SIMD_BLOCK; l++)
            ratio[l] = ratiom[l] = SET_ZERO_R;
        CHARM(lc) lc;
        CHARM(shs_lc_init)(&lc);
        REAL_SIMD zeros[SIMD_BLOCK];
        for (l = 0; l < SIMD_BLOCK; l++)
            zeros[l] = SET_ZERO_R;
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


                ratio[l]  = DIV_R(rref, pnt_r[l]);
                ratiom[l] = ratio[l];


                /* Prepare arrays for sectorial Legendre functions */
                CHARM(leg_func_prepare)(uv, ps + l * SIMD_SIZE * nmax,
                                        ips + l * SIMD_SIZE * nmax, dm, nmax);
            }


            for (size_t p = 0; p < SHS_MAX_NPAR * SIMD_BLOCK; p++)
                fi[p] = SET_ZERO_R;


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Apply polar optimization if asked to do so */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, &u[0],
                                                         SIMD_BLOCK, pt))
                    goto UPDATE_RATIOS;


                /* "anm" and "bnm" coefficients for Legendre recurrence
                 * relations and their derivatives */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);
                if (dorder > 0)
                    CHARM(leg_func_enm)(nmax, m, r, ri, enm);


                /* Computation of the lumped coefficients */
                /* --------------------------------------------------------- */
                /* Summation over harmonic degrees.  Note that the symmetry of
                 * Legendre functions cannot be utilized with scattered points,
                 * hence "SET_ZERO_R".  The "lc.a2" and "lc.b2" are not used
                 * with scattered points. */
#undef KERNEL_IO_PARS
#define KERNEL_IO_PARS (nmax, m, shcs, r_eq_rref, anm, bnm, enm,              \
                        &t[0], &u[0], ps, ips, &ratio[0], &zeros[0],          \
                        &ratiom[0], &zeros[0], &zeros[0], dorder, &lc)
                if ((dr == 0) && (dlat == 0) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr0_dlat0_dlon0)
                        KERNEL_IO_PARS;
                }
                else if ((dr == 1) && (dlat == 0) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr1_dlat0_dlon0)
                        KERNEL_IO_PARS;
                }
                else if ((dr == 2) && (dlat == 0) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr2_dlat0_dlon0)
                        KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 1) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr0_dlat1_dlon0)
                        KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 2) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr0_dlat2_dlon0)
                        KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 0) && (dlon == 1))
                {
                    CHARM(shs_point_kernel_dr0_dlat0_dlon1)
                        KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 0) && (dlon == 2))
                {
                    CHARM(shs_point_kernel_dr0_dlat0_dlon2)
                        KERNEL_IO_PARS;
                }
                else if ((dr == 1) && (dlat == 1) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr1_dlat1_dlon0)
                        KERNEL_IO_PARS;
                }
                else if ((dr == 1) && (dlat == 0) && (dlon == 1))
                {
                    CHARM(shs_point_kernel_dr1_dlat0_dlon1)
                        KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 1) && (dlon == 1))
                {
                    CHARM(shs_point_kernel_dr0_dlat1_dlon1)
                        KERNEL_IO_PARS;
                }
                else if ((dr == GRAD_1) && (dlat == GRAD_1) &&
                         (dlon == GRAD_1))
                {
                    CHARM(shs_point_kernel_grad1) KERNEL_IO_PARS;
                }
                else if ((dr == GRAD_2) && (dlat == GRAD_2) &&
                         (dlon == GRAD_2))
                {
                    CHARM(shs_point_kernel_grad2) KERNEL_IO_PARS;
                }
                /* --------------------------------------------------------- */


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


                    idx = l;
                    LC_CS( );


                    if (grad > 0)
                    {
                        idx += SIMD_BLOCK;
                        LC_CS(r);


                        idx += SIMD_BLOCK;
                        LC_CS(p);
                    }


                    if (grad > 1)
                    {
                        idx += SIMD_BLOCK;
                        LC_CS(rr);


                        idx += SIMD_BLOCK;
                        LC_CS(rp);


                        idx += SIMD_BLOCK;
                        LC_CS(pp);
                    }
                }
                /* ......................................................... */


UPDATE_RATIOS:
                for (l = 0; l < SIMD_BLOCK; l++)
                    ratiom[l] = MUL_R(ratiom[l], ratio[l]);


            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            /* Final part of the synthesis */
            for (size_t p = 0; p < npar; p++)
                CHARM(shs_sctr_mulc)(i, npnt, pnt->type, mur, tmp, tmpv,
                                     &fi[p * SIMD_BLOCK], f[p]);


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
        free(anm);
        free(bnm);
        free(enm);
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
