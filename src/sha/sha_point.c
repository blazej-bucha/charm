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
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../leg/leg_func_xnum.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#if CHARM_PARALLEL
#   include <omp.h>
#endif
/* ------------------------------------------------------------------------- */






void CHARM(sha_point)(const CHARM(crd) *pnt, const REAL *f, unsigned long nmax,
                      CHARM(shc) *shcs, CHARM(err) *err)
{
    /* Some trivial initial error checks */
    /* --------------------------------------------------------------------- */
    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Maximum harmonic degree of the analysis (\"nmax\") "
                       "cannot be larger than "
                       "maximum harmonic degree of spherical harmonic "
                       "coefficients (\"shcs->nmax\").");
        return;
    }


    /* ..................................................................... */
    size_t pnt_nlat = pnt->nlat;
    int pnt_type = pnt->type;
    unsigned long nmax_grd;


    /* Get the maximum degree for which the grid in "pnt" was created.  The
     * value is derived from the number of latitudes "pnt_nlat". */
    if (pnt_type == CHARM_CRD_POINTS_GRID_GL)
    {
        /* In case of the Gauss--Legendre grid, it holds that "nmax_grd
         * = pnt_nlat - 1". */
        nmax_grd = pnt_nlat - 1;
    }
    else if ((pnt_type == CHARM_CRD_POINTS_GRID_DH1) ||
             (pnt_type == CHARM_CRD_POINTS_GRID_DH2))
    {
        /* In case of the Driscoll--Healy grid, it holds that "nmax_grd
         * = (pnt_nlat - 2) / 2". */
        nmax_grd = (pnt_nlat - 2) / 2;


        /* Technically, the Driscoll--Healy grids do not meet our conditions
         * for symmetric grids, because the north pole does not have its
         * negative counterpart (the south pole).  However, we know that except
         * for the north pole, the Driscoll--Healy grids *are* symmetric, so
         * that the symmetry property of Legendre functions could be used.  To
         * make our algorithm work as needed, an easy solution is to increase
         * now "pnt_nlat" by "1". */
        pnt_nlat += 1;
    }
    else
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported \"pnt->type\" for spherical "
                       "harmonic analysis of point data values.");
        return;
    }


    /* Now we check whether the maximum harmonic degree "nmax_grd" is large
     * enough to recover harmonics up to degree "nmax" entered by the user. */
    if (nmax > nmax_grd)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The input data grid \"pnt\" was "
                       "created for a maximum degree that "
                       "is not high enough to recover harmonics up "
                       "to the specified maximum degree of the "
                       "analysis.");
        return;
    }


    /* Get the radius of the sphere "r0", on which the analysis will be
     * performed */
    REAL r0 = pnt->r[0];
    for (size_t i = 1; i < pnt->nlat; i++)
    /* In the line above, we intentionally used "pnt->nlat" instead of
     * "pnt_nlat", since the latter may have been modified, depending on the
     * grid type. */
    {
        if (!CHARM(misc_is_nearly_equal)(pnt->r[i], r0, CHARM(glob_threshold)))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "All spherical radii in \"pnt->r\" must be"
                           "equal.");
            return;
        }
    }
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    /* Check whether the number of latitudes is even or odd */
    /* --------------------------------------------------------------------- */
    _Bool even;
    size_t nlatdo;
    if ((pnt_nlat % 2) == 0)
    {
        even = 1; /* The number of latitudes is an even number. The grid does
                   * not contain the zero latitude */
        nlatdo = pnt_nlat / 2;
    }
    else
    {
        even = 0; /* The number of latitudes is an odd number. The grid could
                   * possibly contain the zero latitude */
        nlatdo = (pnt_nlat + 1) / 2;
    }
    /* --------------------------------------------------------------------- */






    /* Get the number of longitudes in "pnt" */
    /* --------------------------------------------------------------------- */
    size_t pnt_nlon = pnt->nlon;
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    int FAILURE_glob = 0;
    REAL *r                 = NULL;
    REAL *ri                = NULL;
    REAL *dm                = NULL;
    REAL *ftmp_in           = NULL;
    FFTW(complex) *ftmp_out = NULL;
    FFTW(plan) plan         = NULL;
    /* --------------------------------------------------------------------- */






    /* Initializations for recurrence relations to compute Legendre functions
     * */
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
        goto FAILURE; }
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






    /* Auxiliary variable entering the computation of the lumped coefficients
     * */
    /* --------------------------------------------------------------------- */
    REAL c = ADDP(0.0);


    /* Here, we have to use the maximum degree "nmax_grd" for which the input
     * grid in "pnt" was created. */
    if (pnt_type == CHARM_CRD_POINTS_GRID_GL)
        c = PI / (REAL)(nmax_grd + 1);
    else if (pnt_type == CHARM_CRD_POINTS_GRID_DH1)
        c = PI / (REAL)(nmax_grd + 1);
    else if (pnt_type == CHARM_CRD_POINTS_GRID_DH2)
        c = PI / (REAL)(2 * nmax_grd + 2);
    /* --------------------------------------------------------------------- */






    /* Initialize all elements of the output coefficients to zero */
    /* --------------------------------------------------------------------- */
    unsigned long nm_count = ((nmax + 2) * (nmax + 1)) / 2;
    memset(shcs->c[0], 0, nm_count * sizeof(REAL));
    memset(shcs->s[0], 0, nm_count * sizeof(REAL));
    /* --------------------------------------------------------------------- */






    /* Create a plan for FFT */
    /* --------------------------------------------------------------------- */
#if CHARM_PARALLEL && FFTW3_OMP
    if (FFTW(init_threads)() == 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFFTWINIT,
                       "FFTW failed to initialize threads.");
        return;
    }


    FFTW(plan_with_nthreads)(omp_get_max_threads());
#endif
    ftmp_in  = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
    if (ftmp_in == NULL)
    {
        FFTW(free)(ftmp_in);
        FAILURE_glob = 1;
        goto FAILURE;
    }
    ftmp_out = (FFTW(complex) *)FFTW(malloc)((pnt_nlon / 2 + 1) *
                                             sizeof(FFTW(complex)));
    if (ftmp_out == NULL)
    {
        FFTW(free)(ftmp_in);
        FFTW(free)(ftmp_out);
        FAILURE_glob = 1;
        goto FAILURE;
    }


    plan = FFTW(plan_dft_r2c_1d)(pnt_nlon, ftmp_in, ftmp_out, FFTW_ESTIMATE);
    if (plan == NULL)
    {
        FFTW(free)(ftmp_in);
        FFTW(free)(ftmp_out);
        FAILURE_glob = 1;
        goto FAILURE;
    }


    FFTW(free)(ftmp_in); FFTW(free)(ftmp_out);
    /* --------------------------------------------------------------------- */






    /* Loop over latitudes */
    /* --------------------------------------------------------------------- */
    size_t imax = nlatdo + even - 1;
    {
        REAL pnm0, pnm1, pnm2;
        int ix, iy, iz, ixy;
        REAL t, u, x, y, z;


        /* ................................................................. */
        int   *ips  = NULL;
        REAL *ps  = NULL;
#if !(CHARM_PARALLEL)
        REAL *anm = NULL;
        REAL *bnm = NULL;
#endif
        REAL *a  = NULL;
        REAL *b  = NULL;
        REAL *a2 = NULL;
        REAL *b2 = NULL;
        REAL *ftmp_in    = NULL;
        FFTW(complex) *ftmp_out = NULL;


        ips =    (int *)calloc(nmax + 1, sizeof(int));
        if (ips == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ps  = (REAL *)calloc(nmax, sizeof(REAL));
        if (ps == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
#if !(CHARM_PARALLEL)
        anm = (REAL *)calloc(nmax + 1, sizeof(REAL));
        if (anm == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        bnm = (REAL *)calloc(nmax + 1, sizeof(REAL));
        if (bnm == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
#endif
        a  = (REAL *)calloc(pnt_nlon / 2 + 1, sizeof(REAL));
        if (a == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b  = (REAL *)calloc(pnt_nlon / 2 + 1, sizeof(REAL));
        if (b == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        a2 = (REAL *)calloc(pnt_nlon / 2 + 1, sizeof(REAL));
        if (a2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b2 = (REAL *)calloc(pnt_nlon / 2 + 1, sizeof(REAL));
        if (b2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ftmp_in = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
        if (ftmp_in == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ftmp_out = (FFTW(complex) *)FFTW(malloc)((pnt_nlon / 2 + 1) *
                                                 sizeof(FFTW(complex)));
        if (ftmp_out == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        /* ................................................................. */



        REAL amp, amm, bmp, bmm;
        REAL cw;
        REAL wlf;
        _Bool symmi;
        _Bool npm_even; /* True if "n + m" is even */


        for (size_t i = 0; i < nlatdo; i++)
        {
            /* Check whether the symmetry property of LFs can be applied */
            /* ------------------------------------------------------------- */
            if ((pnt_type == CHARM_CRD_POINTS_GRID_DH1 ||
                 pnt_type == CHARM_CRD_POINTS_GRID_DH2) && i == 0)
                /* Do not apply the symmetry property at the north pole of the
                 * Driscoll--Healy grids, as they do not have their negative
                 * counterpart (the south pole). */
                symmi = 0;
            else
            {
                if (i < imax)
                    symmi = 1;
                else
                    symmi = 0;
            }
            /* ------------------------------------------------------------- */


            /* Lumped coefficients for the southern hemisphere (including the
             * equator) */
            /* ------------------------------------------------------------- */
            memcpy(ftmp_in, f + i * pnt_nlon, pnt_nlon * sizeof(REAL));
            FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


            cw = c * pnt->w[i];
            for (size_t j = 0; j < (pnt_nlon / 2 + 1); j++)
            {
                a[j]  =  cw * ftmp_out[j][0];
                b[j]  = -cw * ftmp_out[j][1];
            }
            /* ------------------------------------------------------------- */


            /* Lumped coefficients for the northern hemisphere */
            /* ------------------------------------------------------------- */
            if (symmi)
            {
                memcpy(ftmp_in, f + (pnt_nlat - i - 1) * pnt_nlon,
                       pnt_nlon * sizeof(REAL));
                FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


                cw = c * pnt->w[pnt_nlat - i - 1];
                for (size_t j = 0; j < (pnt_nlon / 2 + 1); j++)
                {
                    a2[j]  =  cw * ftmp_out[j][0];
                    b2[j]  = -cw * ftmp_out[j][1];
                }
            }
            /* ------------------------------------------------------------- */


            /* Prepare arrays for sectorial Legendre functions */
            /* ------------------------------------------------------------- */
            t  = SIN(pnt->lat[i]);
            u  = COS(pnt->lat[i]);


            CHARM(leg_func_prepare)(u, ps, ips, dm, nmax);
            /* ------------------------------------------------------------- */


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(nmax, symmi, r, ri, a, b, a2, b2, shcs, t, ps, ips) \
shared(FAILURE_glob, err) \
private(amp, amm, bmp, bmm) \
private(x, ix, y, iy, wlf, ixy, z, iz) \
private(pnm0, pnm1, pnm2, npm_even)
            {
            /* ............................................................. */
            /* An indicator for failed memory initializations on each thread,
             * a private variable. */
            int FAILURE_priv = 0;


            REAL *anm = NULL;
            REAL *bnm = NULL;


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
#pragma omp critical
            {
                /* At this point, "FAILURE_priv" is equal to 1 if any
                 * allocation failed on a particular thread.  Now, we can add
                 * "FAILURE_priv" from all the threads to get the total
                 * number threads on which an allocation failure occurred (the
                 * "FAILURE_glob" variable).  Note that this code block is
                 * executed with one thread at a time only, so
                 * "FAILURE_glob" can safely be overwritten.
                 * */
                FAILURE_glob += FAILURE_priv;
            }


                /* Now we have to wait until all the threads get here. */
#pragma omp barrier
            /* OK, now let's check on each thread whether there is at least one
             * failed memory allocation among the threads. */
            if (FAILURE_glob > 0)
            {
                /* Ooops, there was indeed a memory allocation failure.  So let
                 * the master thread write down the error to the "err"
                 * variable. */
#pragma omp master
                {
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EMEM, CHARM_ERR_MALLOC_FAILURE);
                }


                /* OK, and now all threads go to the "FAILURE_2_parallel"
                 * label to deallocate all the memory that might be allocated
                 * before the allocation failure. */
                goto FAILURE_2_parallel;
            }
            /* ............................................................. */


#pragma omp for schedule(dynamic)
#endif
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Some useful substitutions */
                /* --------------------------------------------------------- */
                if (symmi)
                {
                    amp = a[m] + a2[m];
                    amm = a[m] - a2[m];
                    bmp = b[m] + b2[m];
                    bmm = b[m] - b2[m];
                }
                else
                {
                    amp = a[m];
                    amm = a[m];
                    bmp = b[m];
                    bmm = b[m];
                }
                /* --------------------------------------------------------- */


                /* Computation of spherical harmonic coefficients */
                if (m == 0)
                {
                    /* Zonal harmonics */
                    /* ----------------------------------------------------- */

                    /* P00 */
                    pnm0 = ADDP(1.0);
                    shcs->c[0][0] += pnm0 * amp; /* C00 */


                    /* P10 */
                    if (nmax >= 1)
                    {
                        pnm1 = ROOT3 * t;
                        shcs->c[0][1] += pnm1 * amm; /* C10 */
                    }


                    /* P20, P30, ..., Pnmax,0 */
                    if (nmax >= 2)
                    {
                        for (unsigned long n = 2; n <= nmax; n++)
                        {
                            pnm2 = anm[n] * t * pnm1 - bnm[n] * pnm0;
                            /* C20, C30, ..., Cnmax,0 */
                            if ((n % 2) == 0)
                                shcs->c[0][n] += pnm2 * amp;
                            else
                                shcs->c[0][n] += pnm2 * amm;


                            pnm0 = pnm1;
                            pnm1 = pnm2;
                        }
                    }
                    /* ----------------------------------------------------- */
                }
                else /* Non-zonal harmonics */
                {

                    /* Sectorial harmonics */
                    /* ----------------------------------------------------- */
                    PNM_SECTORIAL_XNUM(x, ix, ps[m - 1], ips[m - 1], pnm0);


                    shcs->c[m][0] += pnm0 * amp; /* Cmm */
                    shcs->s[m][0] += pnm0 * bmp; /* Smm */
                    /* ----------------------------------------------------- */


                    /* Tesseral harmonics */
                    /* ----------------------------------------------------- */
                    if (m < nmax)
                    {
                        PNM_SEMISECTORIAL_XNUM(x, y, ix, iy, wlf, t,
                                               anm[m + 1], pnm1);


                        shcs->c[m][1] += pnm1 * amm; /* Cm+1,m */
                        shcs->s[m][1] += pnm1 * bmm; /* Sm+1,m */


                        /* Loop over degrees */
                        /* ------------------------------------------------- */
                        /* Is "n + m" even?  Since we start the loop with "n
                         * = m + 2", then the parity of the first "m + 2 + m"
                         * is always even.  Then, it changes with every loop
                         * iteration. */
                        npm_even = 1;
                        for (unsigned long n = (m + 2); n <= nmax;
                             n++, npm_even = !npm_even)
                        {
                            /* Compute tesseral Legendre function */
                            PNM_TESSERAL_XNUM(x, y, z,
                                              ix, iy, iz,
                                              ixy, wlf, t,
                                              anm[n], bnm[n], pnm2, continue);


                            /* Cm+2,m, Cm+3,m, ..., Cnmax,m and Sm+2,m, Sm+3,m,
                             * ..., Snmax,m */
                            if (npm_even)
                            {
                                shcs->c[m][n - m] += pnm2 * amp;
                                shcs->s[m][n - m] += pnm2 * bmp;
                            }
                            else
                            {
                                shcs->c[m][n - m] += pnm2 * amm;
                                shcs->s[m][n - m] += pnm2 * bmm;
                            }


                        } /* End of the loop over harmonic degrees */
                        /* ------------------------------------------------- */


                    } /* End of computation of tesseral harmonics */
                    /* ----------------------------------------------------- */


                } /* End of computation of spherical harmonic coefficients */
                /* --------------------------------------------------------- */


            } /* End of the loop over harmonic orders */


#if CHARM_PARALLEL
FAILURE_2_parallel:
            free(anm); free(bnm);
            }
#endif
            /* ------------------------------------------------------------- */


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


FAILURE_1:
        if (FAILURE_glob != 0)
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);


        FFTW(free)(ftmp_in);
        FFTW(free)(ftmp_out);
        free(ips); free(ps);
        free(a); free(b);
        free(a2); free(b2);
#if !(CHARM_PARALLEL)
        free(anm); free(bnm);
#endif


    }
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if (FAILURE_glob != 0)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(r); free(ri); free(dm);
    FFTW(destroy_plan)(plan);
#if CHARM_PARALLEL && FFTW3_OMP
    FFTW(cleanup_threads)();
#else
    FFTW(cleanup)();
#endif
    /* --------------------------------------------------------------------- */






    /* Multiplication of the "shcs->c" and "shcs->s" arrays by the factor of "1
     * / (4 * pi)", normalization and rescaling */
    /* --------------------------------------------------------------------- */
    /* Normalize the coefficients */
    /* ..................................................................... */
    REAL c2 = ADDP(1.0) / (ADDP(4.0) * PI) * (r0 / shcs->mu);
#if CHARM_PARALLEL
#pragma omp parallel for default(none) shared(shcs, nmax, c2)
#endif
    for (unsigned long m = 0; m <= nmax; m++)
    {
        for (unsigned long n = m; n <= nmax; n++)
        {
            shcs->c[m][n - m] *= c2;
            shcs->s[m][n - m] *= c2;
        }
    }
    /* ..................................................................... */


    /* Rescale the coefficients */
    /* ..................................................................... */
    if (!CHARM(misc_is_nearly_equal)(shcs->r, r0, CHARM(glob_threshold)))
    {
        /* This is the desired radius to scale the output coefficients */
        REAL rtmp = shcs->r;


        /* Now set the "shcs->r" to "r0", since the analysis was conducted on
         * a sphere with the radius "r0", which is different from "shcs->r"
         * (otherwise, the program would not enter this "if" condition). */
        shcs->r = r0;


        /* And now rescale the coefficients from "r0", for which the analysis
         * was done, to the desired "rtmp" */
        CHARM(shc_rescale)(shcs, shcs->mu, rtmp, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    return;
}
