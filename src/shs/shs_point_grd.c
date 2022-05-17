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
#include <fftw3.h>
#include "../prec.h"
#include "shs_point_kernel.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "shs_rpows.h"
#include "shs_check_symmi.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_point_grd)(const CHARM(crd) *pnt, const CHARM(shc) *shcs,
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
    if ((pnt_type == CHARM_CRD_POINTS_GRID_DH1) ||
        (pnt_type == CHARM_CRD_POINTS_GRID_DH2))
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
    else if ((pnt_type == CHARM_CRD_POINTS_GRID_DH1) ||
             (pnt_type == CHARM_CRD_POINTS_GRID_DH2) ||
             (pnt_type == CHARM_CRD_POINTS_GRID_GL))
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
        err_tmp = CHARM(misc_arr_chck_symm)(pnt->lat, pnt_nlat, ADDP(0.0),
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


    /* If true, FFT is applied along the latitude parallels */
    _Bool use_fft;


    /* If "pnt" is a user-defined grid, we have to check whether the
     * longitudinal step is constant.  For quadrature grids, the longitudinal
     * step is constant by definition, so no check is needed. */
    if (pnt->type == CHARM_CRD_POINTS_GRID)
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
    REAL dlon = (pnt_nlon > 1) ? pnt->lon[1] - pnt->lon[0] : ADDP(0.0);


    /* Auxiliary constant to be used only in case the PSLR algorithm is
     * applied.  To suppress, a compiler warning, the value is initialized to
     * zero. */
    REAL lon0 = ADDP(0.0);


    /* Length of the lumped coefficients arrays in case FFT will be applied */
    size_t nlc = pnt_nlon / 2 + 1;


    /* Let's check whether FFT can be employed.  Four conditions must be
     * satisfied to allow FFT: i) the number of grid longitudes must be large
     * enough when compared with "nmax", ii) the longitudinal step must be
     * constant, iii) "pnt->lon[0]" must be zero, and iv) "pnt->lon[pnt_nlon
     * - 1] + dlon" must be equal to "2.0 * PI".
     *
     * The second condition has already been checked, so let's do the remaining
     * checks. */
    if ((pnt_nlon - 1) / 2 >= nmax &&
        CHARM(misc_is_nearly_equal)(pnt->lon[0], ADDP(0.0),
                                    CHARM(glob_threshold)) &&
        CHARM(misc_is_nearly_equal)(pnt->lon[pnt_nlon - 1] + dlon,
                                    ADDP(2.0) * PI,
                                    CHARM(glob_threshold)))
        /* Great, FFT can be applied for this grid! */
        use_fft = 1;
    else
        /* Oh no, FFT cannot be applied for this grid. */
        use_fft = 0;


    if (!use_fft)
    {
        /* OK, so no FFT and the longitudinal step is constant */


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
        lc_tmp  = (FFTW(complex) *)FFTW(malloc)(nlc * sizeof(FFTW(complex)));
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


#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(f, shcs, nmax, pnt, pnt_nlat, pnt_nlon, dm, r, ri, nlatdo, pnt_type) \
shared(lon0, dlon, even, symm, FAILURE_glob, mur, err, nlc, plan, use_fft)
#endif
    {
        /* ................................................................. */
        /* An indicator for failed memory initializations on each thread,
         * a private variable. */
        int FAILURE_priv = 0;


        int  *ips          = NULL;
        REAL *ps           = NULL;
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


        if (use_fft)
        {
            ftmp = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
            if (ftmp == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }


            lc = (FFTW(complex) *)FFTW(malloc)(nlc * sizeof(FFTW(complex)));
            if (lc == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
            memset(lc, 0, nlc * sizeof(FFTW(complex)));
        }
        else
        {
            fi  = (REAL *)calloc(pnt_nlon, sizeof(REAL));
            if (fi == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
        }


        if (symm)
        {
            rpows2 = (REAL *)calloc(nmax + 1, sizeof(REAL));
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


                lc2 = (FFTW(complex) *)FFTW(malloc)(nlc *
                                                    sizeof(FFTW(complex)));
                if (lc2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
                memset(lc2, 0, nlc * sizeof(FFTW(complex)));
            }
            else
            {
                fi2 = (REAL *)calloc(pnt_nlon, sizeof(REAL));
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


        REAL t, u;

        REAL a, b, a2, b2;
        a = b = a2 = b2 = ADDP(0.0);
        REAL c;
        REAL lontmp, clontmp, slontmp, dm0, dm1, dm2, dm02, dm12, dm22;
        REAL cmdlon2, m_fp;

        _Bool symmi;
        size_t row_idx;


#if CHARM_PARALLEL
#pragma omp for
#endif
        for (size_t i = 0; i < nlatdo; i++)
        {
            /* Check whether the symmetry property of LFs can be applied */
            /* ------------------------------------------------------------- */
            symmi = CHARM(shs_check_symmi)(pnt_type, nlatdo, symm, even, i);
            /* ------------------------------------------------------------- */


            /* Prepare arrays for sectorial Legendre functions */
            /* ------------------------------------------------------------- */
            t  = SIN(pnt->lat[i]);
            u  = COS(pnt->lat[i]);


            CHARM(leg_func_prepare)(u, ps, ips, dm, nmax);
            /* ------------------------------------------------------------- */


            /* Pre-compute the powers of "shcs->r / pnt->r[i]" */
            /* ------------------------------------------------------------- */
            CHARM(shs_rpows)(shcs->r, pnt->r[i], rpows, nmax);


            if (symmi)
                CHARM(shs_rpows)(shcs->r, pnt->r[pnt_nlat - i - 1], rpows2,
                                 nmax);
            /* ------------------------------------------------------------- */


            if (use_fft)
            {
                /* Reset the lumped coefficients.  Required in some cases. */
                /* --------------------------------------------------------- */
                memset(lc, 0, nlc * sizeof(FFTW(complex)));


                if (symmi)
                    memset(lc2, 0, nlc * sizeof(FFTW(complex)));
                /* --------------------------------------------------------- */
            }
            else
            {
                /* The "fi" vector represents the synthesized quantity "f" for
                 * the "i"th latitude parallel. Therefore, it needs to be
                 * reinitialized to zero for each "ith" latitude. The same
                 * holds true for "fi2" in case of symmetric grids. */
                /* --------------------------------------------------------- */
                memset(fi, 0, pnt_nlon * sizeof(REAL));


                if (symmi)
                    memset(fi2, 0, pnt_nlon * sizeof(REAL));
                /* --------------------------------------------------------- */
            }


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Computation of the lumped coefficients */
                CHARM(shs_point_kernel)(nmax, m, shcs, anm, bnm,
                                        t, ps, ips, rpows, rpows2,
                                        symmi,
                                        &a, &b, &a2, &b2);


                /* --------------------------------------------------------- */
                if (use_fft)
                {
                    /* Let's prepare the complex Fourier coefficients */
                    c = (m == 0) ? ADDP(1.0) : ADDP(0.5);
                    lc[m][0] =  a * c;
                    lc[m][1] = -b * c;


                    if (symmi)
                    {
                        lc2[m][0] =  a2 * c;
                        lc2[m][1] = -b2 * c;
                    }
                }
                else
                {
                    /* The PSLR algorithm from Balmino et al. (2012) */
                    /* ----------------------------------------------------- */
                    /* The first longitude from the "pnt->lon" array */
                    /* ..................................................... */
                    m_fp    = (REAL)m;
                     lontmp = m_fp * lon0;
                    clontmp = COS(lontmp);
                    slontmp = SIN(lontmp);
                    dm0     = a * clontmp + b * slontmp;
                    fi[0]  += dm0;


                    if (symmi)
                    {
                        dm02    = a2 * clontmp + b2 * slontmp;
                        fi2[0] += dm02;
                    }


                    if (pnt_nlon == 1)
                        continue;
                    /* ..................................................... */


                    /* The second longitude from the "pnt->lon" array */
                    /* ..................................................... */
                     lontmp = m_fp * (lon0 + dlon);
                    clontmp = COS(lontmp);
                    slontmp = SIN(lontmp);
                    dm1     = a * clontmp + b * slontmp;
                    fi[1]  += dm1;


                    if (symmi)
                    {
                        dm12    = a2 * clontmp + b2 * slontmp;
                        fi2[1] += dm12;
                    }


                    if (pnt_nlon == 2)
                        continue;
                    /* ..................................................... */


                    /* The third and all the remaining longitudes from the
                     * "pnt->lon" array */
                    /* ..................................................... */
                    cmdlon2 = ADDP(2.0) * COS(m_fp * dlon);
                    for (size_t j = 2; j < pnt_nlon; j++)
                    {
                        dm2    = cmdlon2 * dm1 - dm0;
                        fi[j] += dm2;
                        dm0    = dm1;
                        dm1    = dm2;
                    }


                    if (symmi)
                    {
                        for (size_t j = 2; j < pnt_nlon; j++)
                        {
                            dm22    = cmdlon2 * dm12 - dm02;
                            fi2[j] += dm22;
                            dm02    = dm12;
                            dm12    = dm22;
                        }
                    }
                    /* ..................................................... */
                    /* ----------------------------------------------------- */
                }


            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            if (use_fft)
            {
                /* Fourier transform along the latitude parallels */
                /* --------------------------------------------------------- */
                FFTW(execute_dft_c2r)(plan, lc, ftmp);
                for (size_t j = 0; j < pnt_nlon; j++)
                    f[i * pnt_nlon + j] = mur * ftmp[j];


                if (symmi)
                {
                    FFTW(execute_dft_c2r)(plan, lc2, ftmp2);
                    for (size_t j = 0; j < pnt_nlon; j++)
                        f[(pnt_nlat - i - 1) * pnt_nlon + j] = mur * ftmp2[j];
                }
                /* --------------------------------------------------------- */
            }
            else
            {
                /* Save the synthesized values to the output "f" array */
                /* --------------------------------------------------------- */
                row_idx = i * pnt_nlon;
                for (size_t j = 0; j < pnt_nlon; j++)
                    f[row_idx + j] = mur * fi[j];


                if (symmi) /* The other hemisphere computed using the symmetry
                            * property of Legendre functions with respect to
                            * the equator */
                {
                    row_idx = (pnt_nlat - i - 1) * pnt_nlon;
                    for (size_t j = 0; j < pnt_nlon; j++)
                        f[row_idx + j] = mur * fi2[j];
                }
                /* --------------------------------------------------------- */
            }


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


        /* Free the heap memory */
        /* ----------------------------------------------------------------- */
FAILURE_2_parallel:
        free(ips); free(ps);
        free(anm); free(bnm);
        free(rpows); free(rpows2);
        FFTW(free)(lc); FFTW(free)(lc2);
        FFTW(free)(ftmp); FFTW(free)(ftmp2);
        free(fi); free(fi2);
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
    if (use_fft)
    {
        FFTW(destroy_plan)(plan);
        FFTW(cleanup)();
    }
    /* --------------------------------------------------------------------- */






    return;
}
