/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "../prec.h"
#include "shs_cell_kernel.h"
#include "shs_grd_lr.h"
#include "shs_grd_lr2.h"
#include "shs_grd_cell_fft_check.h"
#include "shs_grd_fft_lc.h"
#include "shs_cell_check_grd_lons.h"
#include "shs_grd_fft.h"
#include "shs_lc_struct.h"
#include "shs_lc_init.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_gm_hm.h"
#include "../leg/leg_pol_en_fn.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../crd/crd_grd_check_symm.h"
#include "../crd/crd_check_cells.h"
#include "../misc/misc_is_nearly_equal.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
#include "shs_cell_grd.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_grd)(const CHARM(cell) *cell, const CHARM(shc) *shcs,
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
    _Bool even = (cell_nlat % 2) == 0;
    /* ..................................................................... */


    /* Determine whether the cells are symmetric with respect to the
     * equator in the latitudinal direction */
    /* ..................................................................... */
    /* If the grid is symmetric with respect to the equator, then "symm = 1",
     * otherwise "symm = 0". If "symm == 1", the function automatically
     * exploits the symmetry property of Legendre functions in order to
     * accelerate the computation. */


     /* If there is only one cell in the latitudinal direction within the grid,
      * the grid is automatically considered as not symmetric.
      *
      * If there is more than one cell in the latitudinal direction in the
      * grid, let's start by assuming that the grid is symmetric with respect
      * to the equator and check whether this is indeed true */
    _Bool symm = cell_nlat > 1;


    for (size_t i = 0; i < cell_nlat; i++)
    {
        if (CHARM(misc_is_nearly_equal)(cell->latmin[i],
                                        -cell->latmax[cell_nlat - i - 1],
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
    if (symm)
        nlatdo = (cell_nlat + 1 - even) / 2;
    else
        nlatdo = cell_nlat;
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    /* Check the longitudes of the computational grid */
    /* --------------------------------------------------------------------- */
    size_t cell_nlon = cell->nlon;


    REAL deltalon;
    CHARM(shs_cell_check_grd_lons)(cell, &deltalon, err);
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
    size_t nfc = cell_nlon / 2 + 1;


    /* If true, FFT is applied along the latitude parallels */
    _Bool use_fft = CHARM(shs_grd_cell_fft_check)(cell, nmax);
    if (!use_fft)
    {
        /* OK, so no FFT and the longitudinal step is constant */


        /* Get the origin of the longitude "lon" array (will be necessary later
         * for the PSLR algorithm).  The "lon" array is not actually stored in
         * memory, but is given as "(cell->lonmin[j] + cell->lonmax[j]) / 2.0"
         * for "j = 0, 1, ..., nlon - 1". */
        lon0 = (cell->lonmin[0] + cell->lonmax[0]) / PREC(2.0);
    }
    /* --------------------------------------------------------------------- */






    /* Check cell boundaries */
    /* --------------------------------------------------------------------- */
    CHARM(crd_check_cells)(cell, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    int FAILURE_glob   = 0;
    REAL *r            = NULL;
    REAL *ri           = NULL;
    REAL *dm           = NULL;
    REAL *en           = NULL;
    REAL *fn           = NULL;
    REAL *gm           = NULL;
    REAL *hm           = NULL;
    FFTWC(complex) *fc = NULL;
    REAL *ftmp         = NULL;
    FFTW(plan) plan    = NULL;
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
        fc  = (FFTWC(complex) *)FFTW(malloc)(nfc * SIMD_SIZE *
                                             sizeof(FFTWC(complex)));
        if (fc == NULL)
        {
            FFTW(free)(fc);
            FAILURE_glob = 1;
            goto FAILURE;
        }
        ftmp = (REAL *)FFTW(malloc)(cell_nlon * sizeof(REAL));
        if (ftmp == NULL)
        {
            FFTW(free)(fc);
            FFTW(free)(ftmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }


        plan = FFTW(plan_dft_c2r_1d)(cell_nlon, fc, ftmp, FFTW_ESTIMATE);
        if (plan == NULL)
        {
            FFTW(free)(fc);
            FFTW(free)(ftmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }


        FFTW(free)(ftmp); FFTW(free)(fc);
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    REAL mur = shcs->mu / shcs->r;


    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


    /* Radius of the reference sphere that is associated with the spherical
     * harmonic coefficients */
    REAL_SIMD rref = SET1_R(shcs->r);


#if CHARM_OPENMP
#pragma omp parallel default(none) \
shared(f, shcs, nmax, cell, cell_nlat, cell_nlon, dm, en, fn, gm, hm, r, ri) \
shared(nlatdo, lon0, deltalon, even, symm, FAILURE_glob, mur, err, cell_type) \
shared(pt, nfc, plan, use_fft, rref)
#endif
    {
        /* ................................................................. */
        /* An indicator for failed memory initializations on each thread,
         * a private variable. */
        int FAILURE_priv = 0;


        size_t nfi = cell_nlon * SIMD_SIZE * SIMD_BLOCK_S;


        INT *ips1           = NULL;
        INT *ips2           = NULL;
        REAL *ps1           = NULL;
        REAL *ps2           = NULL;
        REAL *latminv       = NULL;
        REAL *latmaxv       = NULL;
        REAL *t1v           = NULL;
        REAL *t2v           = NULL;
        REAL *u1v           = NULL;
        REAL *u2v           = NULL;
        REAL *symmv         = NULL;
        REAL *latsinv       = NULL;
        REAL *fc_simd       = NULL;
        REAL *fc2_simd      = NULL;
        REAL *anm           = NULL;
        REAL *bnm           = NULL;
        REAL *fi            = NULL;
        REAL *fi2           = NULL;
        REAL *ftmp          = NULL;
        FFTWC(complex) *fc  = NULL;
        FFTWC(complex) *fc2 = NULL;
        REAL *cell_rv       = NULL;
        REAL *cell_r2v      = NULL;


        ips1 = (INT *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                            nmax * SIMD_SIZE, sizeof(INT));
        if (ips1 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        ips2 = (INT *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                            nmax * SIMD_SIZE, sizeof(INT));
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
            fc = (FFTWC(complex) *)FFTW(malloc)(nfc * SIMD_SIZE *
                                                sizeof(FFTWC(complex)));
            if (fc == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
            memset(fc, 0, nfc * SIMD_SIZE * sizeof(FFTWC(complex)));
            fc_simd = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                    nfc * SIMD_SIZE * 2,
                                                    sizeof(REAL));
            if (fc_simd == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
        }
        else
        {
            fi = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nfi,
                                               sizeof(REAL));
            if (fi == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
        }


        cell_rv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        if (cell_rv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        if (symm)
        {
            if (use_fft)
            {
                fc2 = (FFTWC(complex) *)FFTW(malloc)(nfc * SIMD_SIZE *
                                                     sizeof(FFTWC(complex)));
                if (fc2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
                memset(fc2, 0, nfc * SIMD_SIZE * sizeof(FFTWC(complex)));
                fc2_simd = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                         nfc * SIMD_SIZE * 2,
                                                         sizeof(REAL));
                if (fc2_simd == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
            }
            else
            {
                fi2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nfi,
                                                    sizeof(REAL));
                if (fi2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
            }


            cell_r2v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                     sizeof(REAL));
            if (cell_r2v == NULL)
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


        REAL_SIMD latmin, latmax;
        REAL_SIMD t1, t2, u1, u2;
        REAL_SIMD symm_simd;


        /* Here, we have pairs of radii-related variables (e.g., "ratio" and
         * "ratio2") not because of the cell boundaries (constant "cell->r" is
         * assumed for each cell), but because of the possible grid symmetry */
        REAL_SIMD cell_r, cell_r2;
        REAL_SIMD ratio, ratio2, ratiom, ratio2m;
        ratio = ratio2 = ratiom = ratio2m = SET_ZERO_R;


        CHARM(lc) lc;
        CHARM(shs_lc_init)(&lc);
        REAL_SIMD a, b, a2, b2;
        a = b = a2 = b2 = SET_ZERO_R;
        REAL_SIMD imm0, imm1, imm2;


        size_t ipv;


        size_t i;
#if CHARM_OPENMP
#pragma omp for schedule(dynamic) private(i)
#endif
        for (i = 0; i < SIMD_MULTIPLE(nlatdo, SIMD_SIZE);
             i += SIMD_SIZE)
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
                    latmaxv[v] = cell->latmax[ipv];
                    latminv[v] = cell->latmin[ipv];
                    t1v[v]     = SIN(latminv[v]);
                    u1v[v]     = COS(latminv[v]);
                    t2v[v]     = SIN(latmaxv[v]);
                    u2v[v]     = COS(latmaxv[v]);
                    cell_rv[v] = cell->r[ipv];


                    if (symm)
                        cell_r2v[v] = cell->r[cell_nlat - ipv - 1];
                }
                else
                {
                    latminv[v] = latmaxv[v] = t1v[v] = u1v[v] = t2v[v] =
                        u2v[v] = cell_rv[v] = PREC(0.0);
                    if (symm)
                        cell_r2v[v] = PREC(0.0);


                    continue;
                }
                /* --------------------------------------------------------- */
            }


            t1        = LOAD_R(&t1v[0]);
            t2        = LOAD_R(&t2v[0]);
            u1        = LOAD_R(&u1v[0]);
            u2        = LOAD_R(&u2v[0]);
            latmin    = LOAD_R(&latminv[0]);
            latmax    = LOAD_R(&latmaxv[0]);
            cell_r    = LOAD_R(&cell_rv[0]);
            symm_simd = LOAD_R(&symmv[0]);


            ratio  = DIV_R(rref, cell_r);
            ratiom = ratio;
            if (symm)
            {
                cell_r2 = LOAD_R(&cell_r2v[0]);
                ratio2  = DIV_R(rref, cell_r2);
                ratio2m = ratio2;
            }


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
                 * have to therefore reset "fc" and "fc2" to zeros.*/
                memset(fc, 0, nfc * SIMD_SIZE * sizeof(FFTWC(complex)));


                if (symm)
                    memset(fc2, 0, nfc * SIMD_SIZE * sizeof(FFTWC(complex)));
            }
            else
            {
                /* The "fi" vector represents the synthesized quantity "f" for
                 * the "ipv"th latitude parallel. Therefore, it needs to be
                 * reinitialized to zero for each "ith" latitude. The same
                 * holds true for "fi2" in case of symmetric grids. */
                memset(fi, 0, nfi * sizeof(REAL));


                if (symm)
                    memset(fi2, 0, nfi * sizeof(REAL));
            }
            /* ------------------------------------------------------------- */


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Apply polar optimization if asked to do so.  Since "u1
                 * = sin(latmin)" and "u2 = sin(latmax)", it is sufficient to
                 * check "u1" only.  In other words, if the polar optimization
                 * can be applied for "u1", it can surely be applied for
                 * "u2". */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, &u1, 1, pt))
                    goto UPDATE_RATIOS;


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
                                       en, fn, gm, hm, ri,
                                       ratio, ratio2,
                                       ratiom, ratio2m,
                                       symm_simd,
                                       &a, &b, &a2, &b2);


                /* The two function calls that follow require "CHARM(lc)" as an
                 * input, so it is prepared here.  We did not want to use
                 * "CHARM(lc)", however, in the previous function call.  The
                 * reason is that "CHARM(lc)" assumes "SIMD_BLOCK_S" larger
                 * than "1" but with cells, the block size is always "1".  In
                 * that case, it may be suboptimal to create "CHARM(lc)" for
                 * some "BLOCK_SIZE > 1", but use it only for a single block
                 * (lots of useless memory jumps that may reduce cache
                 * efficiency).  In the two function calls that follow, the
                 * memory jumps are, however, not at all critical, so they can
                 * rely on "CHARM(lc)" without deteriorating the
                 * performance. */
                lc.a[0]  = a;
                lc.b[0]  = b;
                lc.a2[0] = a2;
                lc.b2[0] = b2;


                if (use_fft)
                    CHARM(shs_grd_fft_lc)(m, deltalon, 0, &lc,
                                          symm, &symm_simd, cell_type,
                                          nfc, fc_simd, fc2_simd);
                else
                    CHARM(shs_grd_lr)(m, lon0, deltalon, cell_nlon, cell_type,
                                      0, nfi, &lc, symm, fi, fi2);


UPDATE_RATIOS:
                ratiom = MUL_R(ratiom, ratio);
                if (symm)
                    ratio2m = MUL_R(ratio2m, ratio2);


            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            if (use_fft)
                /* Fourier transform along the latitude parallels */
                CHARM(shs_grd_fft)(i, cell_type, cell_nlat, cell_nlon,
                                   latsinv, latminv, latmaxv, deltalon,
                                   fc, fc2, nfc, fc_simd, fc2_simd,
                                   mur, plan, symmv, ftmp, f);
            else
                CHARM(shs_grd_lr2)(i, latsinv,
                                   cell_type, cell_nlat, cell_nlon,
                                   symmv, mur, latminv, latmaxv, deltalon,
                                   fi, fi2, f);


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


FAILURE_2_parallel:
        CHARM(free_aligned)(ips1);     CHARM(free_aligned)(ps1);
        CHARM(free_aligned)(ips2);     CHARM(free_aligned)(ps2);
        CHARM(free_aligned)(t1v);      CHARM(free_aligned)(t2v);
        CHARM(free_aligned)(u1v);      CHARM(free_aligned)(u2v);
        CHARM(free_aligned)(cell_rv);  CHARM(free_aligned)(cell_r2v);
        CHARM(free_aligned)(symmv);    CHARM(free_aligned)(latsinv);
        CHARM(free_aligned)(latminv);  CHARM(free_aligned)(latmaxv);
        CHARM(free_aligned)(fi);       CHARM(free_aligned)(fi2);
        free(anm);                     free(bnm);
        FFTW(free)(fc);                FFTW(free)(fc2);
        FFTW(free)(ftmp);
        CHARM(free_aligned)(fc_simd);  CHARM(free_aligned)(fc2_simd);


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
