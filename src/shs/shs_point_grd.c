/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "../prec.h"
#include "shs_rpows.h"
#include "shs_grd_fft.h"
#include "shs_grd_point_fft_check.h"
#include "shs_grd_fft_lc.h"
#include "shs_grd_lr.h"
#include "shs_grd_lr2.h"
#include "shs_point_kernel.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../crd/crd_grd_check_symm.h"
#include "../misc/misc_arr_chck_lin_incr.h"
#include "../misc/misc_arr_chck_symm.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_point_grd)(const CHARM(point) *pnt, const CHARM(shc) *shcs,
                          unsigned long nmax, REAL *f, CHARM(err) *err)
{
    /* Check the latitudes */
    /* --------------------------------------------------------------------- */

    /* Check whether the number of latitudes is even or odd */
    /* ..................................................................... */
    int pnt_type = pnt->type;
    size_t pnt_nlat = pnt->nlat;


    /* If the number of latitudes is even, then "even = 1", otherwise "even =
     * 0" */
    _Bool even;
    if ((pnt_type == CHARM_CRD_POINT_GRID_DH1) ||
        (pnt_type == CHARM_CRD_POINT_GRID_DH2))
    {
        /* For the Driscoll--Healy grids, "pnt_nlat" is always an even number
         * and the grid is not symmetric in terms of our definition (the north
         * pole does not have its negative counterpart -- the south pole).
         * However, we know that except for the north pole, the Driscoll--Healy
         * grids *are* symmetric, so the symmetry property of Legendre
         * functions could be used if the north pole is treated properly.  To
         * this end, let's increase the number of points in "pnt_nlat", so that
         * we can apply our algorithm for symmetric grids.  After increasing
         * "pnt_nlat", its value is odd, so set "even" to zero. */
        pnt_nlat += 1;
        even = 0;
    }
    else
        /* The Gauss--Legendre grid or a user-defined grid */
        even = ((pnt_nlat % 2) == 0) ? 1 : 0;
    /* ..................................................................... */


    /* Determine whether the latitudes are symmetric with respect to the
     * equator */
    /* ..................................................................... */
    /* If the grid is symmetric with respect to the equator, then "symm = 1",
     * otherwise "symm = 0".  If "symm == 1", the function automatically
     * exploits the symmetry property of Legendre functions in order to
     * accelerate the computation. */
    _Bool symm;
    int err_tmp;
    if (pnt_nlat == 1)
        /* The grid is automatically considered as non-symmetric if there is
         * only a single latitude */
        symm = 0;
    else if ((pnt_type == CHARM_CRD_POINT_GRID_DH1) ||
             (pnt_type == CHARM_CRD_POINT_GRID_DH2) ||
             (pnt_type == CHARM_CRD_POINT_GRID_GL))
    {
        /* The Driscoll--Healy grids or the Gauss--Legendre grid.
         *
         * The Gauss--Legendre grid is symmetric by definition, so no check is
         * needed.
         *
         * The Driscoll--Healy grids do not meet our conditions for grids to be
         * symmetric (the south pole is missing).  However, we know that if we
         * omit their north pole, they become symmetric, too.  So let's
         * consider the Driscoll--Healy grids as symmetric by treating the
         * north pole separately later. */
        symm = 1;
    }
    else
    {
        /* User-defined grid, let's do the check */
        err_tmp = CHARM(misc_arr_chck_symm)(pnt->lat, pnt_nlat, PREC(0.0),
                                            CHARM(glob_threshold2), err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }


        symm = (err_tmp == 0) ? 1 : 0;
    }
    /* ..................................................................... */


    /* Finally, if the grid is symmetric, we modify the number of latitudes
     * "pnt_nlat" to be equal to the number of latitudes on one hemisphere
     * only (including the equator if present). This is because the
     * time-consuming "for loop" over evaluation points now needs to run for
     * one hemisphere only, while the results for the other hemisphere are
     * obtained by exploiting the symmetry property of Legendre functions. This
     * reduces the number of Legendre functions that need to be evaluated by
     * a factor of ~2, so saves some computational time */
    /* ..................................................................... */
    size_t nlatdo;
    if (symm)
        nlatdo = (pnt_nlat + 1 - even) / 2;
    else
        nlatdo = pnt_nlat;
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    /* Check the longitudes.  If possible, FFT is employed along the
     * latitudinal parallels.  Otherwise, the PSLR algorithm is used.  The
     * latter is slower, but can be used for any grid with a constant
     * longitudinal sampling.  Below, we determined whether FFT can be applied
     * or not. */
    /* --------------------------------------------------------------------- */
    size_t pnt_nlon = pnt->nlon;


    /* If "pnt" is a user-defined grid with more than one longitude, we have to
     * check whether the longitudinal step is constant.  For quadrature grids,
     * the longitudinal step is constant by definition, so no check is
     * needed. */
    if ((pnt->type == CHARM_CRD_POINT_GRID) && (pnt_nlon > 1))
    {
        err_tmp = CHARM(misc_arr_chck_lin_incr)(pnt->lon, pnt_nlon,
                                                0, 1, CHARM(glob_threshold2),
                                                err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }


        if (err_tmp != 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__,  __func__,
                           CHARM_EFUNCARG,
                           "\"pnt->lon\" is not a linearly increasing "
                           "array within the \"threshold2\".");
            return;
        }
    }


    /* Get the longitudinal step of the grid.  At this point, we know the
     * longitudinal step is constant over each latitude parallel, so let's get
     * the step from, say, the first two longitudes. */
    REAL dlon = (pnt_nlon > 1) ? pnt->lon[1] - pnt->lon[0] : PREC(0.0);


    /* Auxiliary constant to be used only in case the PSLR algorithm is
     * applied.  To suppress, a compiler warning, the value is initialized to
     * zero. */
    REAL lon0 = PREC(0.0);


    /* Length of the lumped coefficients arrays in case FFT will be applied */
    size_t nlc = pnt_nlon / 2 + 1;


    _Bool use_fft = CHARM(shs_grd_point_fft_check)(pnt, dlon, nmax);
    if (!use_fft)
    {
        /* Get the origin of the longitude "pnt->lon" vector (will be necessary
         * later for the PSLR algorithm) */
        lon0 = pnt->lon[0];
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    int FAILURE_glob      = 0;
    REAL *r               = NULL;
    REAL *ri              = NULL;
    REAL *dm              = NULL;
    FFTW(complex) *lc_tmp = NULL;
    REAL *ftmp            = NULL;
    FFTW(plan) plan       = NULL;
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






    /* Create a FFT plan */
    /* --------------------------------------------------------------------- */
    if (use_fft)
    {
        lc_tmp = (FFTW(complex) *)FFTW(malloc)(nlc * sizeof(FFTW(complex)));
        if (lc_tmp == NULL)
        {
            FFTW(free)(lc_tmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }
        ftmp = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
        if (ftmp == NULL)
        {
            FFTW(free)(lc_tmp);
            FFTW(free)(ftmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }


        plan = FFTW(plan_dft_c2r_1d)(pnt_nlon, lc_tmp, ftmp, FFTW_ESTIMATE);
        if (plan == NULL)
        {
            FFTW(free)(lc_tmp);
            FFTW(free)(ftmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }


        FFTW(free)(lc_tmp); FFTW(free)(ftmp);
    }
    /* --------------------------------------------------------------------- */






    /* Loop over grid latitudes */
    /* --------------------------------------------------------------------- */
    REAL mur = shcs->mu / shcs->r;


    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(f, shcs, nmax, pnt, pnt_nlat, pnt_nlon, dm, r, ri, nlatdo, pnt_type) \
shared(lon0, dlon, even, symm, FAILURE_glob, mur, err, nlc, plan, use_fft, pt)
#endif
    {
        /* ................................................................. */
        /* An indicator for failed memory initializations on each thread,
         * a private variable. */
        int FAILURE_priv = 0;


        int  *ips          = NULL;
        REAL *ps           = NULL;
        REAL *latv         = NULL;
        REAL *tv           = NULL;
        REAL *uv           = NULL;
        REAL *symmv        = NULL;
        REAL *latsinv      = NULL;
        REAL *lc_simd      = NULL;
        REAL *lc2_simd     = NULL;
        REAL *anm          = NULL;
        REAL *bnm          = NULL;
        REAL *rpows        = NULL;
        REAL *rpows2       = NULL;
        FFTW(complex) *lc  = NULL;
        FFTW(complex) *lc2 = NULL;
        REAL *ftmp         = NULL;
        REAL *ftmp2        = NULL;
        REAL *fi           = NULL;
        REAL *fi2          = NULL;


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
        latv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE, 
                                             sizeof(REAL));
        if (latv == NULL)
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
        symmv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE, 
                                              sizeof(REAL));
        if (symmv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        latsinv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE, 
                                                sizeof(REAL));
        if (latsinv == NULL)
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


        if (use_fft)
        {
            ftmp = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
            if (ftmp == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
            lc = (FFTW(complex) *)FFTW(malloc)(nlc * SIMD_SIZE *
                                               sizeof(FFTW(complex)));
            if (lc == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
            memset(lc, 0, nlc * SIMD_SIZE * sizeof(FFTW(complex)));
            lc_simd = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                    nlc * SIMD_SIZE * 2,
                                                    sizeof(REAL));
            if (lc_simd == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
        }
        else
        {
            fi = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                               pnt_nlon * SIMD_SIZE,
                                               sizeof(REAL));
            if (fi == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
        }


        if (symm)
        {
            rpows2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                   (nmax + 1) * SIMD_SIZE,
                                                   sizeof(REAL));
            if (rpows2 == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }


            if (use_fft)
            {
                ftmp2 = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
                if (ftmp2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
                lc2 = (FFTW(complex) *)FFTW(malloc)(nlc * SIMD_SIZE *
                                                    sizeof(FFTW(complex)));
                if (lc2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
                memset(lc2, 0, nlc * SIMD_SIZE * sizeof(FFTW(complex)));
                lc2_simd = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                         nlc * SIMD_SIZE * 2,
                                                         sizeof(REAL));
                if (lc2_simd == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
            }
            else
            {
                fi2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                    pnt_nlon * SIMD_SIZE,
                                                    sizeof(REAL));
                if (fi2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
            }
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
            if (CHARM(err_isempty)(err))
                CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                               CHARM_ERR_MALLOC_FAILURE);


            /* OK, and now all threads go to the "FAILURE_2_parallel" label
             * to deallocate all the memory that might be allocated before the
             * allocation failure. */
            goto FAILURE_2_parallel;
        }
        /* ................................................................. */


        REAL_SIMD t, u, symm_simd;
        REAL_SIMD a, b, a2, b2;
        a = b = a2 = b2 = SET_ZERO_R;

        size_t ipv;


#if CHARM_PARALLEL
#pragma omp for schedule(dynamic)
#endif
        for (size_t i = 0; i < SIMD_GET_MULTIPLE(nlatdo); i += SIMD_SIZE)
        {
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                /* Check whether the symmetry property of LFs can be applied */
                /* --------------------------------------------------------- */
                ipv = i + v;
                CHARM(crd_grd_check_symm)(ipv, v, pnt_type, nlatdo, symm, even,
                                          symmv, latsinv);


                if (latsinv[v] == 1)
                {
                    latv[v] = pnt->lat[ipv];
                    tv[v]   = SIN(latv[v]);
                    uv[v]   = COS(latv[v]);
                }
                else
                {
                    latv[v] = tv[v] = uv[v] = PREC(0.0);
                    continue;
                }
                /* --------------------------------------------------------- */


                /* Pre-compute the powers of "shcs->r / pnt->r[ipv]" */
                /* --------------------------------------------------------- */
                CHARM(shs_rpows)(v, shcs->r, pnt->r[ipv], rpows, nmax);


                if (symmv[v])
                    CHARM(shs_rpows)(v, shcs->r, pnt->r[pnt_nlat - ipv - 1],
                                     rpows2, nmax);
                /* --------------------------------------------------------- */
            }


            t         = LOAD_R(&tv[0]);
            u         = LOAD_R(&uv[0]);
            symm_simd = LOAD_R(&symmv[0]);


            /* Prepare arrays for sectorial Legendre functions */
            /* ------------------------------------------------------------- */
            CHARM(leg_func_prepare)(uv, ps, ips, dm, nmax);
            /* ------------------------------------------------------------- */


            if (use_fft)
            {
                /* Reset the lumped coefficients.  Required in some cases. */
                /* --------------------------------------------------------- */
                memset(lc, 0, nlc * SIMD_SIZE * sizeof(FFTW(complex)));


                if (symm)
                    memset(lc2, 0, nlc * SIMD_SIZE * sizeof(FFTW(complex)));
                /* --------------------------------------------------------- */
            }
            else
            {
                /* The "fi" vector represents the synthesized quantity "f" for
                 * the "ipv"th latitude parallel. Therefore, it needs to be
                 * reinitialized to zero for each "ith" latitude. The same
                 * holds true for "fi2" in case of symmetric grids. */
                /* --------------------------------------------------------- */
                memset(fi, 0, pnt_nlon * SIMD_SIZE * sizeof(REAL));


                if (symm)
                    memset(fi2, 0, pnt_nlon * SIMD_SIZE * sizeof(REAL));
                /* --------------------------------------------------------- */
            }


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


                /* Computation of the lumped coefficients */
                CHARM(shs_point_kernel)(nmax, m, shcs, anm, bnm,
                                        t, ps, ips, rpows, rpows2,
                                        symm_simd, &a, &b, &a2, &b2);


                /* --------------------------------------------------------- */
                if (use_fft)
                    CHARM(shs_grd_fft_lc)(m, dlon, a, b, a2, b2, symm,
                                          symm_simd, pnt_type,
                                          lc_simd, lc2_simd);
                else
                    CHARM(shs_grd_lr)(m, lon0, dlon, pnt_nlon, pnt_type,
                                      a, b, a2, b2, symm, fi, fi2);
                /* --------------------------------------------------------- */


            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            if (use_fft)
            {
                /* Fourier transform along the latitude parallels */
                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    if (latsinv[v] == 0)
                        continue;


                    CHARM(shs_grd_fft)(i, v, pnt_nlat, pnt_nlon, lc, lc2, nlc,
                                       lc_simd, lc2_simd, mur, plan, symmv,
                                       ftmp, ftmp2, f);
                }
            }
            else
            {
                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    if (latsinv[v] == 0)
                        continue;


                    CHARM(shs_grd_lr2)(i, v, pnt_nlat, pnt_nlon, symmv, mur,
                                       fi, fi2, f);
                }
                /* --------------------------------------------------------- */
            }


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


        /* Free the heap memory */
        /* ----------------------------------------------------------------- */
FAILURE_2_parallel:
        CHARM(free_aligned)(ips);       CHARM(free_aligned)(ps);
        CHARM(free_aligned)(latv);
        CHARM(free_aligned)(tv);        CHARM(free_aligned)(uv);
        CHARM(free_aligned)(symmv);     CHARM(free_aligned)(latsinv);
        CHARM(free_aligned)(rpows);     CHARM(free_aligned)(rpows2);
        CHARM(free_aligned)(fi);        CHARM(free_aligned)(fi2);
        CHARM(free_aligned)(lc_simd);   CHARM(free_aligned)(lc2_simd);
        free(anm); free(bnm);
        FFTW(free)(lc); FFTW(free)(lc2);
        FFTW(free)(ftmp); FFTW(free)(ftmp2);
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
    if (use_fft)
    {
        FFTW(destroy_plan)(plan);
        FFTW(cleanup)();
    }
    /* --------------------------------------------------------------------- */






    return;
}
