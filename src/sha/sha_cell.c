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
#include "../shc/shc_reset_coeffs.h"
#include "../shs/shs_cell_check_grd_lons.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_gm_hm.h"
#include "../leg/leg_pol_en_fn.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../leg/leg_func_xnum.h"
#include "../crd/crd_grd_check_symm.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../misc/misc_is_nearly_equal.h"
#if CHARM_PARALLEL
#   include <omp.h>
#endif
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
/* ------------------------------------------------------------------------- */






void CHARM(sha_cell)(const CHARM(crd) *cell, const REAL *f, unsigned long nmax,
                     int method, CHARM(shc) *shcs, CHARM(err) *err)
{
    /* Some trivial initial error checks */
    /* --------------------------------------------------------------------- */
    if (method != CHARM_SHA_CELL_AQ)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported \"method\" of spherical "
                       "harmonic analysis of block-mean data values.");
        return;
    }


    if (cell->type != CHARM_CRD_CELLS_GRID)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported \"cell->type\" for spherical "
                       "harmonic analysis of block-mean data values.");
        return;
    }


    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Maximum harmonic degree of the analysis "
                       "(\"nmax\") cannot be larger than "
                       "maximum harmonic degree of spherical harmonic "
                       "coefficients (\"shcs->nmax\").");
        return;
    }


    /* Get the radius of the sphere "r0", on which the analysis will be
     * performed */
    REAL r0 = cell->r[0];
    size_t cell_nlat = cell->nlat;
    for (size_t i = 1; i < cell_nlat; i++)
    {
        if (!CHARM(misc_is_nearly_equal)(cell->r[i], r0,
                                         CHARM(glob_threshold)))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "All spherical radii in \"cell->r\" must be"
                           "equal.");
            return;
        }
    }
    /* --------------------------------------------------------------------- */






    /* Check the latitudes of the computational grid */
    /* --------------------------------------------------------------------- */
    /* The first latitude must be equal to "-PI_2". */
    /* ..................................................................... */
    if (CHARM(misc_is_nearly_equal)(cell->lat[0], -PI_2,
                                    CHARM(glob_threshold)) != 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lat[0]\" has to be equal to \"-PI / 2.0\" "
                       "within the \"threshold\".");
        return;
    }
    /* ..................................................................... */


    /* Check whether we have enough data in the latitudinal direction */
    /* ..................................................................... */
    if (cell_nlat < nmax + 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Not enough data in the latitudinal direction.  "
                       "The required minimum is \"nmax + 1\".");
        return;
    }
    /* ..................................................................... */


    /* The last latitude must be equal to "PI_2". */
    /* ..................................................................... */
    if (CHARM(misc_is_nearly_equal)(cell->lat[2 * cell_nlat - 1],
                                    PI_2, CHARM(glob_threshold)) != 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lat[2 * cell->nlat - 1]\" has to "
                       "be equal to \"PI / 2.0\" within the "
                       "\"threshold\".");
        return;
    }
    /* ..................................................................... */


    /* Now check whether the grid cells are nicely aligned in the latitudinal
     * direction. */
    /* ..................................................................... */
    for (size_t i = 0; i < (cell_nlat - 1); i++)
    {
        if (CHARM(misc_is_nearly_equal)(cell->lat[2 * i + 1],
                                        cell->lat[2 * (i + 1)],
                                        CHARM(glob_threshold)) != 1)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The grid cells are not nicely aligned in the "
                           "latitudinal direction.  It must hold that "
                           "\"cell->lat[2 * i + 1] == cell->lat[2 * (i + 1)]\""
                           " for all \"i = 0, 1, ..., cell->nlat - 1\".");
            return;
        }
    }
    /* ..................................................................... */


    /* Check whether the number of cells in the latitudinal direction is even
     * or odd */
    /* ..................................................................... */
    _Bool even;
    if ((cell_nlat % 2) == 0)
        /* The number of cells in the latitudinal direction is an even number
         * */
        even = 1;
    else
        /* The number of cells in the latitudinal direction is an odd number */
        even = 0;
    /* ..................................................................... */


    /* Determine whether the latitudes are symmetric with respect to the
     * equator */
    /* ..................................................................... */
    /* If the grid is symmetric with respect to the equator, then "symm = 1",
     * otherwise "symm = 0". If "symm == 1", the function automatically
     * exploits the symmetry property of Legendre functions in order to
     * accelerate the computation */
    _Bool symm;
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
            /* The grid is not symmetric */
            symm = 0;


            /* Exiting the loop with the final decision: the grid is not
             * symmetric with respect to the equator ("symm = 0") */
            break;
        }
    }
    /* ..................................................................... */


    /* Finally, if the grid is symmetric, we modify the number of cells in the
     * latitudinal direction "cell_nlat" to be equal to the number of cells on
     * one hemisphere only (including the equator if present). This is because
     * the time-consuming "for loop" over evaluation points now needs to run
     * for one hemisphere only, while the results for the other hemisphere are
     * obtained by exploiting the symmetry property of Legendre functions. This
     * reduces the number of Legendre functions that need to be evaluated by
     * a factor of ~2, so saves some computational time */
    /* ..................................................................... */
    size_t nlatdo;
    if (symm)
    {
        if (even == 1)
            nlatdo = cell_nlat / 2;
        else
            nlatdo = (cell_nlat + 1) / 2;
    }
    else
        nlatdo = cell_nlat;
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    /* Check the longitudes of the computational grid */
    /* --------------------------------------------------------------------- */
    /* Check whether we have enough data in the longitudinal direction */
    /* ..................................................................... */
    size_t cell_nlon = cell->nlon;
    if (cell_nlon < 2 * nmax + 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Not enough data in the longitudinal direction.  "
                       "The required minimum is \"2 * nmax + 1\".");
        return;
    }
    size_t cell_nlon_fft = cell_nlon / 2 + 1;
    /* ..................................................................... */


    if (CHARM(misc_is_nearly_equal)(cell->lon[0], PREC(0.0),
                                    CHARM(glob_threshold)) != 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lon[0]\" has to be equal to \"0.0\" "
                       "within the \"threshold\".");
        return;
    }


    if (CHARM(misc_is_nearly_equal)(cell->lon[2 * cell->nlon - 1],
                                    PREC(2.0) * PI,
                                    CHARM(glob_threshold)) != 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lon[2 * cell->nlon - 1]\" has to be "
                       "equal to \"2.0 * PI\" within the "
                       "\"threshold\".");
        return;
    }


    REAL dlon;
    CHARM(shs_cell_check_grd_lons)(cell, &dlon, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    /* Now check whether the grid cells are nicely aligned in the longitudinal
     * direction.  Due to the previous error checks, this can be done
     * only with the first grid cell in the longitudinal direction. */
    if (cell_nlon > 1)
    {
        if (CHARM(misc_is_nearly_equal)(cell->lon[2], cell->lon[1],
                                        CHARM(glob_threshold)) == 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The grid cells are not nicely aligned in the "
                           "longitudinal direction.  It must hold that "
                           "\"cell->lon[2 * j] == cell->lon[2 * j - 1]\" for "
                           "all \"j = 1, 2, ..., cell->nlon - 1\".");
            return;
        }
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    int FAILURE_glob = 0;
    REAL *r  = NULL;
    REAL *ri = NULL;
    REAL *dm = NULL;
    REAL *en = NULL;
    REAL *fn = NULL;
    REAL *gm = NULL;
    REAL *hm = NULL;
    REAL *ftmp_in           = NULL;
    FFTW(complex) *ftmp_out = NULL;
    FFTW(plan) plan         = NULL;
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






    /* Set all coefficients in "shcs" to zero */
    CHARM(shc_reset_coeffs)(shcs);






    /* Create a plan for FFT */
    /* --------------------------------------------------------------------- */
#if CHARM_PARALLEL && FFTW3_OMP
    int err_fftw = FFTW(init_threads)();
    if (err_fftw == 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__,  __func__, CHARM_EFFTWINIT,
                       "FFTW failed to initialize threads.");
        return;
    }

    FFTW(plan_with_nthreads)(omp_get_max_threads());
#endif
    ftmp_in  = (REAL *)FFTW(malloc)(cell_nlon * sizeof(REAL));
    if (ftmp_in == NULL)
    {
        FFTW(free)(ftmp_in);
        FAILURE_glob = 1;
        goto FAILURE;
    }
    ftmp_out = (FFTW(complex) *)FFTW(malloc)(sizeof(FFTW(complex)) * 
                                             cell_nlon_fft);
    if (ftmp_out == NULL)
    {
        FFTW(free)(ftmp_in);
        FFTW(free)(ftmp_out);
        FAILURE_glob = 1;
        goto FAILURE;
    }


    plan = FFTW(plan_dft_r2c_1d)(cell_nlon, ftmp_in, ftmp_out, FFTW_ESTIMATE);
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
    {
        REAL_SIMD x1, y1, z1, t1, u1;
        REAL_SIMD x2, y2, z2, t2, u2;
        RI_SIMD   ix1, iy1, iz1, ixy1;
        RI_SIMD   ix2, iy2, iz2, ixy2;
        REAL_SIMD pnm0_latmin, pnm1_latmin, pnm2_latmin;
        REAL_SIMD pnm0_latmax, pnm1_latmax, pnm2_latmax;
        REAL_SIMD latmin, latmax;
#ifdef SIMD
        RI_SIMD    zero_ri = SET_ZERO_RI;
        RI_SIMD    one_ri  = SET1_RI(1);
        RI_SIMD    mone_ri = SET1_RI(-1);
        REAL_SIMD  zero_r  = SET_ZERO_R;
        REAL_SIMD  BIG_r   = SET1_R(BIG);
        REAL_SIMD  BIGI_r  = SET1_R(BIGI);
        REAL_SIMD  BIGS_r  = SET1_R(BIGS);
        REAL_SIMD  BIGSI_r = SET1_R(BIGSI);
        REAL_SIMD  tmp1_r,  tmp2_r;
        MASK_SIMD  mask1, mask2;
        MASK2_SIMD mask3;
        ABS_R_INIT;
#endif
        REAL_SIMD symm_simd, latsin;
        REAL_SIMD am, bm, a2m, b2m;
        REAL_SIMD cm_simd, sm_simd;
        REAL_SIMD amp, amm, bmp, bmm;
        REAL_SIMD in0, inm0, inm1, inm2;
        REAL_SIMD w;
        REAL cm, sm, mr;
        _Bool npm_even; /* True if "n + m" is even */
        size_t ipv; /* "i + v" */
        unsigned long nmm; /* "n - m" */


        /* ................................................................. */
        int  *ips1    = NULL;
        int  *ips2    = NULL;
        REAL *ps1     = NULL;
        REAL *ps2     = NULL;
        REAL *latminv = NULL;
        REAL *latmaxv = NULL;
        REAL *t1v     = NULL;
        REAL *t2v     = NULL;
        REAL *u1v     = NULL;
        REAL *u2v     = NULL;
        REAL *symmv   = NULL;  /* To be store only "1.0" or "0.0" */
        REAL *latsinv = NULL;  /* To be store only "1.0" or "0.0" */
        REAL *a       = NULL;
        REAL *b       = NULL;
        REAL *a2      = NULL;
        REAL *b2      = NULL;
#if !(CHARM_PARALLEL)
        REAL *anm     = NULL;
        REAL *bnm     = NULL;
#endif
        REAL *imm     = NULL;
        REAL *ftmp_in = NULL;
        FFTW(complex) *ftmp_out = NULL;


        ips1 = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE * nmax, 
                                            sizeof(int));
        if (ips1 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ips2 = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE * nmax, 
                                            sizeof(int));
        if (ips2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ps1 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE * nmax, 
                                            sizeof(REAL));
        if (ps1 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ps2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE * nmax, 
                                            sizeof(REAL));
        if (ps2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        latminv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        if (latminv == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        latmaxv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        if (latmaxv == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        t1v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (t1v == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        t2v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (t2v == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        u1v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (u1v == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        u2v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
        if (u2v == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        symmv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE, 
                                              sizeof(REAL));
        if (symmv == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        latsinv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE, 
                                                sizeof(REAL));
        if (latsinv == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        a = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, 
                                          SIMD_SIZE * cell_nlon_fft, 
                                          sizeof(REAL));
        if (a == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, 
                                          SIMD_SIZE * cell_nlon_fft, 
                                          sizeof(REAL));
        if (b == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        a2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, 
                                           SIMD_SIZE * cell_nlon_fft, 
                                           sizeof(REAL));
        if (a2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, 
                                           SIMD_SIZE * cell_nlon_fft, 
                                           sizeof(REAL));
        if (b2 == NULL)
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
        imm = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, 
                                            SIMD_SIZE * (nmax + 1), 
                                            sizeof(REAL));
        if (imm == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ftmp_in = (REAL *)FFTW(malloc)(cell_nlon * sizeof(REAL));
        if (ftmp_in == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ftmp_out = (FFTW(complex) *)FFTW(malloc)(sizeof(FFTW(complex)) *
                                                 cell_nlon_fft);
        if (ftmp_out == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        /* ................................................................. */



        for (size_t i = 0; i < SIMD_GET_MULTIPLE(nlatdo); i += SIMD_SIZE)
        {
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                /* Check whether the symmetry property of LFs needs to be
                 * applied */
                /* --------------------------------------------------------- */
                ipv = i + v;
                CHARM(crd_grd_check_symm)(ipv, v, cell->type, nlatdo, symm,
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


                /* Lumped coefficients for the southern hemisphere (including
                 * the equator) */
                /* --------------------------------------------------------- */
                memcpy(ftmp_in, f + ipv * cell_nlon, cell_nlon * sizeof(REAL));
                FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


                for (size_t j = 0; j < cell_nlon_fft; j++)
                {
                    a[j * SIMD_SIZE + v] =  ftmp_out[j][0];
                    b[j * SIMD_SIZE + v] = -ftmp_out[j][1];
                }
                /* --------------------------------------------------------- */


                /* Lumped coefficients for the northern hemisphere */
                /* --------------------------------------------------------- */
                if (symmv[v])
                {
                    memcpy(ftmp_in, f + (cell_nlat - ipv - 1) * cell_nlon,
                           cell_nlon * sizeof(REAL));
                    FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


                    for (size_t j = 0; j < cell_nlon_fft; j++)
                    {
                        a2[j * SIMD_SIZE + v] =  ftmp_out[j][0];
                        b2[j * SIMD_SIZE + v] = -ftmp_out[j][1];
                    }
                }
                /* --------------------------------------------------------- */
            }


            t1        = LOAD_R(&t1v[0]);
            t2        = LOAD_R(&t2v[0]);
            u1        = LOAD_R(&u1v[0]);
            u2        = LOAD_R(&u2v[0]);
            latmin    = LOAD_R(&latminv[0]);
            latmax    = LOAD_R(&latmaxv[0]);
            symm_simd = LOAD_R(&symmv[0]);
            latsin    = LOAD_R(&latsinv[0]);


            /* Prepare arrays for sectorial Legendre functions */
            /* ------------------------------------------------------------- */
            CHARM(leg_func_prepare)(u1v, ps1, ips1, dm, nmax);
            CHARM(leg_func_prepare)(u2v, ps2, ips2, dm, nmax);
            /* ------------------------------------------------------------- */


            /* Pre-compute the sectorial "imm" integrals.  This is necessary
             * for the "defined(CHARM_PARALLEL)" parallelization strategy, but
             * can be used (and in fact it is) also with the other
             * parallelization strategies. */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 1; m <= nmax; m++)
            {
                /* Pmm for "latmin" */
#ifdef SIMD
                PNM_SECTORIAL_XNUM_SIMD(x1, ix1,
                                        ps1[SIMD_SIZE * (m - 1)],
                                        ips1[SIMD_SIZE * (m - 1)],
                                        pnm0_latmin, BIG_r, zero_r, zero_ri,
                                        mone_ri, mask1,  mask2, SECTORIALS1);
#else
                PNM_SECTORIAL_XNUM(x1, ix1, ps1[m - 1], ips1[m - 1],
                                   pnm0_latmin);
#endif


                /* Pmm for "latmax" */
#ifdef SIMD
                PNM_SECTORIAL_XNUM_SIMD(x2, ix2,
                                        ps2[SIMD_SIZE * (m - 1)],
                                        ips2[SIMD_SIZE * (m - 1)],
                                        pnm0_latmax, BIG_r, zero_r, zero_ri,
                                        mone_ri, mask1, mask2, SECTORIALS2);
#else
                PNM_SECTORIAL_XNUM(x2, ix2, ps2[m - 1], ips2[m - 1],
                                   pnm0_latmax);
#endif


                /* Imm */
                /* ..................................................... */
                if (m == 1) /* I11 */
                {
                    STORE_R(&imm[m * SIMD_SIZE],
                            MUL_R(SET1_R(ROOT3 / PREC(2.0)),
                                  SUB_R(SUB_R(MUL_R(t2, u2),
                                              SUB_R(SET1_R(PI_2), latmax)),
                                        SUB_R(MUL_R(t1, u1),
                                              SUB_R(SET1_R(PI_2), latmin)))));
                }
                else if (m == 2) /* I22 */
                {
                    STORE_R(&imm[m * SIMD_SIZE], 
                            MUL_R(SET1_R(SQRT(PREC(15.0)) / PREC(6.0)),
                                  SUB_R(MUL_R(t2, SUB_R(SET1_R(PREC(3.0)),
                                                        MUL_R(t2, t2))),
                                        MUL_R(t1, SUB_R(SET1_R(PREC(3.0)),
                                                        MUL_R(t1, t1))))));
                }
                else /* I33, I44, ..., Inmax,nmax */
                {
                    STORE_R(&imm[m * SIMD_SIZE], 
                            ADD_R(MUL_R(SET1_R(gm[m]),
                                        LOAD_R(&imm[(m - 2) * SIMD_SIZE])),
                                  MUL_R(SET1_R(PREC(1.0) / (REAL)(m + 1)),
                                        SUB_R(MUL_R(t2, pnm0_latmax),
                                              MUL_R(t1, pnm0_latmin)))));
                }
                /* ..................................................... */
            }
            /* ------------------------------------------------------------- */






            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
#if CHARM_PARALLEL
#   ifdef SIMD
        #pragma omp parallel default(none) \
            shared(nmax, t1, t2, u1, u2, symm_simd) \
            shared(shcs, en, fn, gm, hm, imm, ps1, ps2, ips1, ips2, r, ri) \
            shared(a, b, a2, b2, dlon, latsin) \
            shared(FAILURE_glob, err) \
            private(am, bm, a2m, b2m) \
            private(cm, sm, x1, x2, ix1, ix2, y1, y2, iy1, iy2, w, mr) \
            private(z1, z2, iz1, iz2, ixy1, ixy2) \
            private(amp, amm, bmp, bmm, cm_simd, sm_simd) \
            private(pnm0_latmin, pnm0_latmax, pnm1_latmin, pnm1_latmax) \
            private(pnm2_latmin, pnm2_latmax, in0, inm0, inm1, inm2) \
            private(npm_even, nmm) \
            shared(zero_ri, one_ri, mone_ri, zero_r) \
            shared(BIG_r, BIGI_r, BIGS_r, BIGSI_r, ABS_R_MASK) \
            private(tmp1_r, tmp2_r, mask1, mask2, mask3)
#   else
        #pragma omp parallel default(none) \
            shared(nmax, t1, t2, u1, u2, symm_simd) \
            shared(shcs, en, fn, gm, hm, imm, ps1, ps2, ips1, ips2, r, ri) \
            shared(a, b, a2, b2, dlon, latsin) \
            shared(FAILURE_glob, err) \
            private(am, bm, a2m, b2m) \
            private(cm, sm, x1, x2, ix1, ix2, y1, y2, iy1, iy2, w, mr) \
            private(z1, z2, iz1, iz2, ixy1, ixy2) \
            private(amp, amm, bmp, bmm, cm_simd, sm_simd) \
            private(pnm0_latmin, pnm0_latmax, pnm1_latmin, pnm1_latmax) \
            private(pnm2_latmin, pnm2_latmax, in0, inm0, inm1, inm2) \
            private(npm_even, nmm)
#   endif
            {
            /* ............................................................. */
            /* An indicator for failed memory initializations on each
             * thread, a private variable. */
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
                 * allocation failed on a particular thread.  Now, we can
                 * add "FAILURE_priv" from all the threads to get the
                 * total number threads on which an allocation failure
                 * occurred (the "FAILURE_glob" variable).  Note that
                 * this code block is executed with one thread at a time
                 * only, so "FAILURE_glob" can safely be overwritten.
                 * */
                FAILURE_glob += FAILURE_priv;
            }


            /* Now we have to wait until all the threads get here. */
#pragma omp barrier
            /* OK, now let's check on each thread whether there is at least
             * one failed memory allocation among the threads. */
            if (FAILURE_glob > 0)
            {
                /* Ooops, there was indeed a memory allocation failure.  So
                 * let the master thread write down the error to the "err"
                 * variable. */
#pragma omp master
                if (CHARM(err_isempty)(err))
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EMEM, CHARM_ERR_MALLOC_FAILURE);


                /* OK, and now all threads go to the
                 * "FAILURE_2_parallel" label to deallocate all the
                 * memory that might be allocated before the allocation
                 * failure. */
                goto FAILURE_2_parallel;
            }
            /* ............................................................. */


#pragma omp for schedule(dynamic)
#endif
            for (unsigned long m = 0; m <= nmax; m++)
            {


                /* Due to the use of the mean values, some additional terms
                 * need to be taken into account when compared with the
                 * harmonic analysis based on point values (see, e.g., Colombo,
                 * 1981) */
                /* --------------------------------------------------------- */
                am  = LOAD_R(&a[SIMD_SIZE * m]);
                bm  = LOAD_R(&b[SIMD_SIZE * m]);
                a2m = LOAD_R(&a2[SIMD_SIZE * m]);
                b2m = LOAD_R(&b2[SIMD_SIZE * m]);


                if (m == 0)
                {
                    amp = MUL_R(MUL_R(ADD_R(am, MUL_R(symm_simd, a2m)),
                                      SET1_R(dlon)), latsin);
                    amm = MUL_R(MUL_R(SUB_R(am, MUL_R(symm_simd, a2m)),
                                      SET1_R(dlon)), latsin);
                }
                else
                {
                    /* Useful substitution */
                    mr = (REAL)m;
                    cm = (COS(mr * dlon) - PREC(1.0)) / mr;
                    sm = SIN(mr * dlon) / mr;
                    cm_simd = SET1_R(cm);
                    sm_simd = SET1_R(sm);


                    amp = MUL_R(ADD_R(MUL_R(ADD_R(am, MUL_R(symm_simd, a2m)),
                                            sm_simd),
                                      MUL_R(ADD_R(bm, MUL_R(symm_simd, b2m)),
                                            cm_simd)), latsin);
                    amm = MUL_R(ADD_R(MUL_R(SUB_R(am, MUL_R(symm_simd, a2m)),
                                            sm_simd),
                                      MUL_R(SUB_R(bm, MUL_R(symm_simd, b2m)),
                                            cm_simd)), latsin);
                    bmp = MUL_R(SUB_R(MUL_R(ADD_R(bm, MUL_R(symm_simd, b2m)),
                                            sm_simd),
                                      MUL_R(ADD_R(am, MUL_R(symm_simd, a2m)),
                                            cm_simd)), latsin);
                    bmm = MUL_R(SUB_R(MUL_R(SUB_R(bm, MUL_R(symm_simd, b2m)),
                                            sm_simd),
                                      MUL_R(SUB_R(am, MUL_R(symm_simd, a2m)),
                                            cm_simd)), latsin);
                }
                /* --------------------------------------------------------- */



                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Computation of spherical harmonic coefficients */
                if (m == 0)
                {






                    /* Zonal Legendre functions and their integrals */
                    /* ----------------------------------------------------- */
                    /* P00 */
                    /* ..................................................... */
                    pnm0_latmin = SET1_R(PREC(1.0));
                    pnm0_latmax = SET1_R(PREC(1.0));
                    /* ..................................................... */


                    /* P10 */
                    /* ..................................................... */
                    pnm1_latmin = t1;
                    pnm1_latmax = t2;
                    /* ..................................................... */


                    /* I00 */
                    /* ..................................................... */
                    in0 = SUB_R(t2, t1);
                    /* ..................................................... */


                    /* Spherical harmonic coefficients */
                    /* ..................................................... */
                    shcs->c[0][0] += SUM_R(MUL_R(in0, amp));
                    /* ..................................................... */


                    if ((nmax + 1) >= 2)
                    {
                        /* Is "n + m" even?  Since we start the loop with "n
                         * = 1" and "m = 0", then the parity of the first "n
                         * + m" is always odd.  Then, it changes with every
                         * loop iteration. */
                        npm_even = 0;
                        for (unsigned long n = 1; n <= nmax;
                             n++, npm_even = !npm_even)
                        {
                            /* P20, P30, ..., Pnmax+1,0
                             * Since we use locally the "n" variable to compute
                             * Legendre polynomials of degree "n + 1", the "n
                             * + 1"th elements have to be taken from the
                             * vectors "en" and "fn" below. Once again, note
                             * that when "m == 0", then "pnm0", "pnm1" and
                             * "pnm2" for latmini and latmaxi represent
                             * un-normalized Legendre polynomials. These are
                             * needed to get the integrals of fully-normalized
                             * Legendre polynomials */
                            /* ............................................. */
                            pnm2_latmin = SUB_R(MUL_R(SET1_R(en[n + 1]),
                                                      MUL_R(t1, pnm1_latmin)),
                                                MUL_R(SET1_R(fn[n + 1]),
                                                      pnm0_latmin));
                            pnm2_latmax = SUB_R(MUL_R(SET1_R(en[n + 1]),
                                                      MUL_R(t2, pnm1_latmax)),
                                                MUL_R(SET1_R(fn[n + 1]),
                                                      pnm0_latmax));
                            /* ............................................. */


                            /* I10, I20, ..., Inmax,0
                             * Computed from Pn+1,0 and Pn-1,0 */
                            /* ............................................. */
                            in0 = MUL_R(SET1_R(ri[2 * n + 1]),
                                        SUB_R(SUB_R(pnm2_latmax, pnm0_latmax),
                                              SUB_R(pnm2_latmin,
                                                    pnm0_latmin)));
                            /* ............................................. */


                            /* Spherical harmonic coefficients */
                            /* ............................................. */
                            /* C10, C20, ..., Cnmax,0 */
                            if (npm_even)
                                shcs->c[0][n] += SUM_R(MUL_R(in0, amp));
                            else
                                shcs->c[0][n] += SUM_R(MUL_R(in0, amm));
                            /* ............................................. */


                            pnm0_latmin = pnm1_latmin;
                            pnm1_latmin = pnm2_latmin;


                            pnm0_latmax = pnm1_latmax;
                            pnm1_latmax = pnm2_latmax;

                        }
                    }
                    /* ----------------------------------------------------- */







                }
                else /* Non-zonal harmonics */
                {






                    /* Sectorial Legendre functions and their integrals */
                    /* ----------------------------------------------------- */

                    /* Pmm for "latmin" */
#ifdef SIMD
                    PNM_SECTORIAL_XNUM_SIMD(x1, ix1,
                                            ps1[SIMD_SIZE * (m - 1)],
                                            ips1[SIMD_SIZE * (m - 1)],
                                            pnm0_latmin, BIG_r, zero_r,
                                            zero_ri, mone_ri, mask1, mask2,
                                            SECTORIALS3);
#else
                    PNM_SECTORIAL_XNUM(x1, ix1, ps1[m - 1], ips1[m - 1],
                                       pnm0_latmin);
#endif


                    /* Pmm for "latmax" */
#ifdef SIMD
                    PNM_SECTORIAL_XNUM_SIMD(x2, ix2,
                                            ps2[SIMD_SIZE * (m - 1)],
                                            ips2[SIMD_SIZE * (m - 1)],
                                            pnm0_latmax, BIG_r, zero_r,
                                            zero_ri, mone_ri, mask1, mask2,
                                            SECTORIALS4);
#else
                    PNM_SECTORIAL_XNUM(x2, ix2, ps2[m - 1], ips2[m - 1],
                                       pnm0_latmax);
#endif


                    /* Imm */
                    /* ..................................................... */
                    inm0 = LOAD_R(&imm[m * SIMD_SIZE]);
                    /* ..................................................... */


                    /* Spherical harmonic coefficients */
                    /* ..................................................... */
                    shcs->c[m][0] += SUM_R(MUL_R(inm0, amp));
                    shcs->s[m][0] += SUM_R(MUL_R(inm0, bmp));
                    /* ..................................................... */
                    /* ----------------------------------------------------- */






                    /* Tesseral harmonics */
                    /* ----------------------------------------------------- */
                    if (m < nmax)
                    {

                        /* Pm+1,m for "latmin" */
#ifdef SIMD
                        PNM_SEMISECTORIAL_XNUM_SIMD(x1, y1, ix1, iy1, w, t1, 
                                                    anm[m + 1], pnm1_latmin, 
                                                    mask1, mask2, mask3, 
                                                    zero_r, zero_ri, mone_ri, 
                                                    BIG_r, BIGS_r,  BIGI_r,
                                                    SEMISECTORIALS1);
#else
                        PNM_SEMISECTORIAL_XNUM(x1, y1, ix1, iy1, w, t1,
                                               anm[m + 1], pnm1_latmin);
#endif


                        /* Pm+1,m for "latmax" */
#ifdef SIMD
                        PNM_SEMISECTORIAL_XNUM_SIMD(x2, y2, ix2, iy2, w, t2, 
                                                    anm[m + 1], pnm1_latmax, 
                                                    mask1, mask2, mask3, 
                                                    zero_r, zero_ri, mone_ri, 
                                                    BIG_r, BIGS_r,  BIGI_r,
                                                    SEMISECTORIALS2);
#else
                        PNM_SEMISECTORIAL_XNUM(x2, y2, ix2, iy2, w, t2,
                                               anm[m + 1], pnm1_latmax);
#endif


                        /* Im+1,m */
                        /* ................................................. */
                        /* This is not a typo, "pnm0" are indeed required here
                         * */
                        inm1 = -MUL_R(SET1_R(anm[m + 1] / (REAL)(m + 2)),
                                      SUB_R(MUL_R(MUL_R(u2, u2), pnm0_latmax),
                                            MUL_R(MUL_R(u1, u1), 
                                                  pnm0_latmin)));
                        /* ................................................. */


                        /* Spherical harmonic coefficients */
                        /* ................................................. */
                        shcs->c[m][1] += SUM_R(MUL_R(inm1, amm));
                        shcs->s[m][1] += SUM_R(MUL_R(inm1, bmm));
                        /* ................................................. */


                        /* Pm+2,m, Pm+3,m, ..., Pnmax,m and their integrals */
                        /* ................................................. */
                        /* Is "n + m" even?  Since we start the loop with "n
                         * = m + 2", then the parity of the first "m + 2 + m"
                         * is always even.  Then, it changes with every loop
                         * iteration. */
                        npm_even = 1;
                        for (unsigned long n = (m + 2); n <= nmax;
                             n++, npm_even = !npm_even)
                        {

                            /* Pm+2,m, Pm+3,m, ..., Pnmax,m for "latmin" */
#ifdef SIMD
                            PNM_TESSERAL_XNUM_SIMD(x1, y1, z1, ix1, iy1, iz1, 
                                                   ixy1,
                                                   w, t1, anm[n], bnm[n], 
                                                   pnm2_latmin,
                                                   tmp1_r, tmp2_r,
                                                   mask1, mask2, 
                                                   mask3, zero_r,
                                                   zero_ri, one_ri,
                                                   BIG_r, BIGI_r,
                                                   BIGS_r, BIGSI_r,
                                                   TESSERALS1, TESSERALS2);
#else
                            PNM_TESSERAL_XNUM(x1, y1, z1,
                                              ix1, iy1, iz1,
                                              ixy1, w, t1,
                                              anm[n], bnm[n],
                                              pnm2_latmin,
                                              pnm2_latmin = PREC(0.0));
#endif


                            /* Pm+2,m, Pm+3,m, ..., Pnmax,m for "latmax" */
#ifdef SIMD
                            PNM_TESSERAL_XNUM_SIMD(x2, y2, z2, ix2, iy2, iz2, 
                                                   ixy2,
                                                   w, t2, anm[n], bnm[n], 
                                                   pnm2_latmax,
                                                   tmp1_r, tmp2_r,
                                                   mask1, mask2, 
                                                   mask3, zero_r,
                                                   zero_ri, one_ri,
                                                   BIG_r, BIGI_r,
                                                   BIGS_r, BIGSI_r,
                                                   TESSERALS3, TESSERALS4);
#else
                            PNM_TESSERAL_XNUM(x2, y2, z2,
                                              ix2, iy2, iz2,
                                              ixy2, w, t2,
                                              anm[n], bnm[n],
                                              pnm2_latmax,
                                              pnm2_latmax = PREC(0.0));
#endif


                            /* Im+2,m, Im+3,m, ..., Inmax,m */
                            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                            /* This is not a typo, "pnm1" are indeed required
                             * here */
                            inm2 = SUB_R(MUL_R(MUL_R(SET1_R(hm[n]),
                                                     SET1_R(bnm[n])), inm0),
                                         MUL_R(SET1_R(anm[n] / (REAL)(n + 1)),
                                               SUB_R(MUL_R(MUL_R(u2, u2),
                                                           pnm1_latmax),
                                                     MUL_R(MUL_R(u1, u1),
                                                           pnm1_latmin))));


                            inm0 = inm1;
                            inm1 = inm2;
                            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


                            /* Spherical harmonic coefficients */
                            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                            /* Cm+2,m, Cm+3,m, ..., Cnmax,m and Sm+2,m, Sm+3,m,
                             * * ..., Snmax,m */
                            nmm = n - m;
                            if (npm_even)
                            {
                                shcs->c[m][nmm] += SUM_R(MUL_R(inm2, amp));
                                shcs->s[m][nmm] += SUM_R(MUL_R(inm2, bmp));
                            }
                            else
                            {
                                shcs->c[m][nmm] += SUM_R(MUL_R(inm2, amm));
                                shcs->s[m][nmm] += SUM_R(MUL_R(inm2, bmm));
                            }
                            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


                            pnm1_latmin = pnm2_latmin;
                            pnm1_latmax = pnm2_latmax;

                        } /* End of the loop over harmonic degrees */
                        /* ------------------------------------------------- */


                    } /* End of computation of tesseral harmonics */
                    /* ----------------------------------------------------- */

                } /* End of computation of spherical harmonic coefficients */
                /* --------------------------------------------------------- */


            } /* End of the loop over harmonic orders */
#if CHARM_PARALLEL
FAILURE_2_parallel:
            free(anm);  free(bnm);
            }
#endif
            /* ------------------------------------------------------------- */


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


FAILURE_1:
        if ((FAILURE_glob != 0) && CHARM(err_isempty)(err))
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);


        FFTW(free)(ftmp_in);          FFTW(free)(ftmp_out);
        CHARM(free_aligned)(t1v);     CHARM(free_aligned)(t2v);
        CHARM(free_aligned)(u1v);     CHARM(free_aligned)(u2v);
        CHARM(free_aligned)(symmv);   CHARM(free_aligned)(latsinv);
        CHARM(free_aligned)(ips1);    CHARM(free_aligned)(ps1);
        CHARM(free_aligned)(ips2);    CHARM(free_aligned)(ps2);
        CHARM(free_aligned)(latminv); CHARM(free_aligned)(latmaxv);
#if !(CHARM_PARALLEL)
        free(anm);  free(bnm);
#endif
        CHARM(free_aligned)(a);    CHARM(free_aligned)(b);
        CHARM(free_aligned)(a2);   CHARM(free_aligned)(b2);
        CHARM(free_aligned)(imm);


    }
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if ((FAILURE_glob != 0) && CHARM(err_isempty)(err))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(r); free(ri); free(dm); free(en); free(fn); free(gm); free(hm);
    FFTW(destroy_plan)(plan);
    FFTW(cleanup)();
#if CHARM_PARALLEL && FFTW3_OMP
    FFTW(cleanup_threads)();
#else
    FFTW(cleanup)();
#endif


    if (FAILURE_glob != 0)
        return;
    /* --------------------------------------------------------------------- */






    /* Multiplication of the "shcs->c" and "shcs->s" arrays by the factor of "1
     * / (4 * pi)", normalization and rescaling */
    /* --------------------------------------------------------------------- */
    /* Normalize the coefficients */
    /* ..................................................................... */
    REAL c2 = PREC(1.0) / (PREC(4.0) * PI) * (r0 / shcs->mu);


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
