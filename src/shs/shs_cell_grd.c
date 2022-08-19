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
#include "shs_cell_kernel.h"
#include "shs_grd_lr.h"
#include "shs_grd_lr2.h"
#include "shs_grd_fft_check.h"
#include "shs_grd_fft_lc.h"
#include "shs_cell_check_grd_lons.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_gm_hm.h"
#include "../leg/leg_pol_en_fn.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "shs_rpows.h"
#include "../crd/crd_grd_check_symm.h"
#include "shs_grd_fft.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_grd)(const CHARM(crd) *cell, const CHARM(shc) *shcs,
                         unsigned long nmax, REAL *f, CHARM(err) *err)
{
    /* Check the latitudes of the computational grid */
    /* --------------------------------------------------------------------- */
    /* Check whether the number of cells in the latitudinal direction in "cell"
     * is even or odd */
    /* ..................................................................... */
    int cell_type = cell->type;
    size_t cell_nlat = cell->nlat;


    /* If the number is even, then "even = 1", otherwise "even = 0" */
    _Bool even = ((cell_nlat % 2) == 0) ? 1 : 0;
    /* ..................................................................... */


    /* Determine whether the cells are symmetric with respect to the
     * equator in the latitudinal direction */
    /* ..................................................................... */
    _Bool symm; /* If the grid is symmetric with respect to the equator, then
                 * "symm = 1", otherwise "symm = 0". If "symm == 1", the
                 * function automatically exploits the symmetry property of
                 * Legendre functions in order to accelerate the computation */
    if (cell_nlat == 1)
        /* If there is only one cell in the latitudinal direction within the
         * grid, the grid is automatically considered as not symmetric */
        symm = 0;
    else
        /* If there is more than one cell in the latitudinal direction in the
         * grid, let's start by assuming that the grid is symmetric with
         * respect to the equator and check whether this is indeed true */
        symm = 1;


    for (size_t i = 0; i < cell_nlat; i++)
    {
        if (CHARM(misc_is_nearly_equal)(cell->lat[i],
                                        -cell->lat[(2 * cell_nlat) - i - 1],
                                        CHARM(glob_threshold2)) == 0)
        {
            symm = 0; /* The grid is not symmetric */

            break; /* Exiting the loop with the final decision: the grid is not
                    * symmetric with respect to the equator ("symm = 0") */
        }
    }
    /* ..................................................................... */


    /* Finally, if the grid is symmetric, we modify the number of latitudes
     * "cell_nlat" to be equal to the number of latitudes on one hemisphere
     * only (including the equator if present). This is because the
     * time-consuming "for loop" over evaluation cells now needs to run for one
     * hemisphere only, while the results for the other hemisphere are obtained
     * by exploiting the symmetry property of Legendre functions. This reduces
     * the number of Legendre functions needed to be evaluated by a factor of
     * ~2, so saves some computational time. */
    /* ..................................................................... */
    size_t nlatdo;
    if (symm == 1)
        nlatdo = (cell_nlat + 1 - even) / 2;
    else
        nlatdo = cell_nlat;
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    /* Check the longitudes of the computational grid */
    /* --------------------------------------------------------------------- */
    size_t cell_nlon = cell->nlon;


    REAL dlon;
    CHARM(shs_cell_check_grd_lons)(cell, &dlon, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    /* Auxiliary constant to be used only in case the PSLR algorithm is
     * applied.  To suppress, a compiler warning, the value is initialized to
     * zero. */
    REAL lon0 = PREC(0.0);


    /* Length of the lumped coefficients arrays in case FFT will be applied */
    size_t nlc = cell_nlon / 2 + 1;


    /* If true, FFT is applied along the latitude parallels */
    _Bool use_fft = CHARM(shs_grd_fft_check)(cell, dlon, nmax);
    if (!use_fft)
    {
        /* OK, so no FFT and the longitudinal step is constant */


        /* Get the origin of the longitude "lon" array (will be necessary later
         * for the PSLR algorithm).  The "lon" array is not actually stored in
         * memory, but is given as "(cell->lon[2 * j] + cell->lon[2 * j + 1])
         * / 2.0" for "j = 0, 1, ..., nlon - 1". */
        lon0 = (cell->lon[0] + cell->lon[1]) / PREC(2.0);
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    int FAILURE_glob  = 0;
    REAL *r           = NULL;
    REAL *ri          = NULL;
    REAL *dm          = NULL;
    REAL *en          = NULL;
    REAL *fn          = NULL;
    REAL *gm          = NULL;
    REAL *hm          = NULL;
    FFTW(complex) *lc = NULL;
    REAL *ftmp        = NULL;
    FFTW(plan) plan   = NULL;
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






    /* Create a FFT plan */
    /* --------------------------------------------------------------------- */
    if (use_fft)
    {
        lc  = (FFTW(complex) *)FFTW(malloc)(nlc * SIMD_SIZE * 
                                            sizeof(FFTW(complex)));
        if (lc == NULL)
        {
            FFTW(free)(lc);
            FAILURE_glob = 1;
            goto FAILURE;
        }
        ftmp = (REAL *)FFTW(malloc)(cell_nlon * sizeof(REAL));
        if (ftmp == NULL)
        {
            FFTW(free)(lc);
            FFTW(free)(ftmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }


        plan = FFTW(plan_dft_c2r_1d)(cell_nlon, lc, ftmp, FFTW_ESTIMATE);
        if (plan == NULL)
        {
            FFTW(free)(lc);
            FFTW(free)(ftmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }


        FFTW(free)(ftmp); FFTW(free)(lc);
    }
    /* --------------------------------------------------------------------- */






    /* Loop over grid latitudes */
    /* --------------------------------------------------------------------- */
    REAL mur = shcs->mu / shcs->r;


#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(f, shcs, nmax, cell, cell_nlat, cell_nlon, dm, en, fn, gm, hm, r, ri) \
shared(nlatdo, lon0, dlon, even, symm, FAILURE_glob, mur, err, cell_type) \
shared(nlc, plan, use_fft)
#endif
    {
        /* ................................................................. */
        /* An indicator for failed memory initializations on each thread,
         * a private variable. */
        int FAILURE_priv = 0;


        int  *ips1         = NULL;
        int  *ips2         = NULL;
        REAL *ps1          = NULL;
        REAL *ps2          = NULL;
        REAL *latminv      = NULL;
        REAL *latmaxv      = NULL;
        REAL *t1v          = NULL;
        REAL *t2v          = NULL;
        REAL *u1v          = NULL;
        REAL *u2v          = NULL;
        REAL *symmv        = NULL;
        REAL *latsinv      = NULL;
        REAL *lc_simd      = NULL;
        REAL *lc2_simd     = NULL;
        REAL *anm          = NULL;
        REAL *bnm          = NULL;
        REAL *fi           = NULL;
        REAL *fi2          = NULL;
        REAL *ftmp         = NULL;
        REAL *ftmp2        = NULL;
        FFTW(complex) *lc  = NULL;
        FFTW(complex) *lc2 = NULL;
        REAL *rpows        = NULL;
        REAL *rpows2       = NULL;


        ips1 = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nmax * SIMD_SIZE,
                                            sizeof(int));
        if (ips1 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        ips2 = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nmax *SIMD_SIZE,
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
        t1v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (t1v == NULL)
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
        u1v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (u1v == NULL)
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


        if (use_fft)
        {
            ftmp = (REAL *)FFTW(malloc)(cell_nlon * sizeof(REAL));
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
                                               cell_nlon * SIMD_SIZE,
                                               sizeof(REAL));
            if (fi == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
        }


        rpows = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                              (nmax + 1) * SIMD_SIZE,
                                              sizeof(REAL));
        if (rpows == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        if (symm)
        {
            if (use_fft)
            {
                ftmp2 = (REAL *)FFTW(malloc)(cell_nlon * sizeof(REAL));
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
                                                    cell_nlon * SIMD_SIZE,
                                                    sizeof(REAL));
                if (fi2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
            }


            rpows2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                   (nmax + 1) * SIMD_SIZE,
                                                   sizeof(REAL));
            if (rpows2 == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
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


        REAL_SIMD latmin, latmax;
        REAL_SIMD t1, t2, u1, u2;
        REAL_SIMD symm_simd;


        REAL_SIMD a, b, a2, b2;
        a = b = a2 = b2 = SET_ZERO_R;
        REAL_SIMD imm0, imm1, imm2;
        REAL dsigma, mur_dsigma;

        size_t ipv;


#if CHARM_PARALLEL
#pragma omp for
#endif
        for (size_t i = 0; i < SIMD_GET_MULTIPLE(nlatdo); i += SIMD_SIZE)
        {
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                /* Check whether the symmetry property of LFs needs to be
                 * exploited */
                /* --------------------------------------------------------- */
                ipv = i + v;
                CHARM(crd_grd_check_symm)(ipv, v, cell_type, nlatdo, symm,
                                          even, symmv, latsinv);


                if (latsinv[v] == 1)
                {
                    latminv[v] = cell->lat[2 * ipv];
                    latmaxv[v] = cell->lat[2 * ipv + 1];
                    t1v[v]     = SIN(latminv[v]);
                    u1v[v]     = COS(latminv[v]);
                    t2v[v]     = SIN(latmaxv[v]);
                    u2v[v]     = COS(latmaxv[v]);
                }
                else
                {
                    latminv[v] = latmaxv[v] = t1v[v] = u1v[v] = t2v[v] = 
                        u2v[v] = PREC(0.0);
                    continue;
                }
                /* --------------------------------------------------------- */


                /* Pre-compute the powers of "shcs->r / cell->r[i + v]" */
                /* --------------------------------------------------------- */
                CHARM(shs_rpows)(v, shcs->r, cell->r[ipv], rpows, nmax);


                if (symmv[v])
                    CHARM(shs_rpows)(v, shcs->r, cell->r[cell_nlat - ipv - 1],
                                     rpows2, nmax);
                /* --------------------------------------------------------- */
            }


            t1        = LOAD_R(&t1v[0]);
            t2        = LOAD_R(&t2v[0]);
            u1        = LOAD_R(&u1v[0]);
            u2        = LOAD_R(&u2v[0]);
            latmin    = LOAD_R(&latminv[0]);
            latmax    = LOAD_R(&latmaxv[0]);
            symm_simd = LOAD_R(&symmv[0]);


            /* Prepare arrays for sectorial Legendre functions */
            /* ------------------------------------------------------------- */
            CHARM(leg_func_prepare)(u1v, ps1, ips1, dm, nmax);
            CHARM(leg_func_prepare)(u2v, ps2, ips2, dm, nmax);
            /* ------------------------------------------------------------- */


            /* ------------------------------------------------------------- */
            if (use_fft)
            {
                /* Important note.  For the "c2r" FFT transform, FFTW
                 * overwrites even the *input* array.  For some combinations of
                 * "nmax" and "cell_nlon", this is not a problem, while for
                 * other combinations, this causes incorrect results.  Here, we
                 * have to therefore reset "lc" and "lc2" to zeros.*/
                memset(lc, 0, nlc * SIMD_SIZE * sizeof(FFTW(complex)));


                if (symm)
                    memset(lc2, 0, nlc * SIMD_SIZE * sizeof(FFTW(complex)));
            }
            else
            {
                /* The "fi" vector represents the synthesized quantity "f" for
                 * the "ipv"th latitude parallel. Therefore, it needs to be
                 * reinitialized to zero for each "ith" latitude. The same
                 * holds true for "fi2" in case of symmetric grids. */
                memset(fi, 0, cell_nlon * SIMD_SIZE * sizeof(REAL));


                if (symm)
                    memset(fi2, 0, cell_nlon * SIMD_SIZE * sizeof(REAL));
            }
            /* ------------------------------------------------------------- */


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {

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
                                       en, fn, gm, hm, ri, rpows, rpows2,
                                       symm_simd,
                                       &a, &b, &a2, &b2);


                if (use_fft)
                    CHARM(shs_grd_fft_lc)(m, dlon, a, b, a2, b2, symm,
                                          symm_simd, cell_type,
                                          lc_simd, lc2_simd);
                else
                    CHARM(shs_grd_lr)(m, lon0, dlon, cell_nlon, cell_type,
                                      a, b, a2, b2, symm, fi, fi2);

            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            if (use_fft)
            {
                /* Fourier transform along the latitude parallels */
                /* --------------------------------------------------------- */
                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    if (latsinv[v] == 0)
                        continue;


                    /* Cell area on the unit sphere and some useful constants
                     * */
                    dsigma     = (SIN(latmaxv[v]) - SIN(latminv[v])) * dlon;
                    mur_dsigma = mur / dsigma;


                    CHARM(shs_grd_fft)(i, v, cell_nlat, cell_nlon, lc, lc2,
                                       nlc, lc_simd, lc2_simd, mur_dsigma,
                                       plan, symmv, ftmp, ftmp2, f);
                }
                /* --------------------------------------------------------- */
            }
            else
            {
                /* Final synthesis */
                /* --------------------------------------------------------- */
                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    if (latsinv[v] == 0)
                        continue;


                    /* Cell area on the unit sphere and some useful constants
                     * */
                    dsigma     = (SIN(latmaxv[v]) - SIN(latminv[v])) * dlon;
                    mur_dsigma = mur / dsigma;


                    CHARM(shs_grd_lr2)(i, v, cell_nlat, cell_nlon, symmv,
                                       mur_dsigma, fi, fi2, f);
                }
                /* --------------------------------------------------------- */
            }


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


FAILURE_2_parallel:
        CHARM(free_aligned)(ips1);     CHARM(free_aligned)(ps1);
        CHARM(free_aligned)(ips2);     CHARM(free_aligned)(ps2);
        CHARM(free_aligned)(t1v);      CHARM(free_aligned)(t2v);
        CHARM(free_aligned)(u1v);      CHARM(free_aligned)(u2v);
        CHARM(free_aligned)(rpows);    CHARM(free_aligned)(rpows2);
        CHARM(free_aligned)(symmv);    CHARM(free_aligned)(latsinv);
        CHARM(free_aligned)(latminv);  CHARM(free_aligned)(latmaxv);
        CHARM(free_aligned)(fi);       CHARM(free_aligned)(fi2);
        free(anm);                     free(bnm);
        FFTW(free)(lc);                FFTW(free)(lc2);
        FFTW(free)(ftmp);              FFTW(free)(ftmp2);
        CHARM(free_aligned)(lc_simd);  CHARM(free_aligned)(lc2_simd);


    } /* End of "#pragma omp parallel" */
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if ((FAILURE_glob != 0) && CHARM(err_isempty)(err))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(r); free(ri);
    free(dm); free(en);
    free(fn); free(gm);
    free(hm);
    if (use_fft)
    {
        FFTW(destroy_plan)(plan);
        FFTW(cleanup)();
    }
    /* --------------------------------------------------------------------- */






    return;
}
