/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <fftw3.h>
#if HAVE_OPENMP
#   include <omp.h>
#endif
#include "../prec.h"
#include "../crd/crd_point_quad_equator.h"
#include "../crd/crd_point_quad_get_nmax_from_nlat.h"
#include "../crd/crd_point_isGLGrid.h"
#include "../crd/crd_point_isDHGrid.h"
#include "../crd/crd_point_isQuadGrid.h"
#include "../shc/shc_reset_coeffs.h"
#include "../shc/shc_block_struct.h"
#include "../shc/shc_block_init.h"
#include "../shc/shc_block_free.h"
#include "../shc/shc_block_set_coeffs.h"
#include "../shc/shc_block_reset_coeffs.h"
#include "../shc/shc_block_get_idx.h"
#include "../shc/shc_block_set_mfirst.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../leg/leg_func_xnum.h"
#include "../leg/leg_func_use_xnum.h"
#include "../crd/crd_grd_check_symm.h"
#include "../crd/crd_point_get_local_nlat.h"
#include "../crd/crd_point_get_local_nlon.h"
#include "../crd/crd_point_get_local_0_start.h"
#include "../shs/shs_get_imax.h"
#if HAVE_MPI
#   include "../mpi/mpi_err_gather.h"
#   include "../mpi/mpi_check_point_shc_err.h"
#endif
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_isempty_all_mpi_processes.h"
#include "../err/err_omp_mpi.h"
#include "../misc/misc_is_nearly_equal.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#include "../misc/misc_sd_calloc.h"
#include "../glob/glob_get_sha_block_lat_multiplier.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
#undef CHECK_NULL
#define CHECK_NULL(x, barrier)                                                \
        if ((x) == NULL)                                                      \
        {                                                                     \
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,     \
                           CHARM_ERR_MALLOC_FAILURE);                         \
            goto barrier;                                                     \
        }






#undef CHECK_NULL_OMP
#define CHECK_NULL_OMP(x, err_priv, barrier)                                  \
        if ((x) == NULL)                                                      \
        {                                                                     \
            err_priv = 1;                                                     \
            goto barrier;                                                     \
        }






#define CS_SUM(cs, pnm, ab)                                                   \
    cs_sum = SET_ZERO_R;                                                      \
    for (l = 0; l < BLOCK_A; l++)                                             \
    {                                                                         \
        cs_sum = ADD_R(cs_sum, MUL_R((pnm)[l], (ab)[l]));                     \
    }                                                                         \
    (cs) += SUM_R(cs_sum);






#define LOOP_ITER(n, a, b)                                                    \
    anms = SET1_R(anm[(n)]);                                                  \
    bnms = SET1_R(bnm[(n)]);                                                  \
                                                                              \
                                                                              \
    for (l = 0; l < BLOCK_A; l++)                                             \
    {                                                                         \
        PNM_RECURRENCE(x[l], y[l], pnm2[l], t[l], anms, bnms);                \
        RECURRENCE_NEXT_ITER(y[l], x[l], pnm2[l]);                            \
    }                                                                         \
                                                                              \
                                                                              \
    CS_SUM(shcs_block->c[idx], pnm2, a);                                      \
    CS_SUM(shcs_block->s[idx++], pnm2, b);
/* ------------------------------------------------------------------------- */






void CHARM(sha_point)(const CHARM(point) *pnt,
                      const REAL *f,
                      unsigned long nmax,
                      CHARM(shc) *shcs,
                      CHARM(err) *err)
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


#if HAVE_MPI
    CHARM(mpi_check_point_shc_err)(pnt, shcs, err);
    if (!CHARM(mpi_err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
#endif


    const int pnt_type = pnt->type;
    if (!CHARM(crd_point_isQuadGrid)(pnt_type))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported \"pnt->type\" for spherical "
                       "harmonic analysis of point data values.");
        return;
    }


    /* ..................................................................... */
    size_t pnt_nlat = CHARM(crd_point_get_local_nlat)(pnt);
    const size_t local_0_start = CHARM(crd_point_get_local_0_start)(pnt);


    /* Get the radius of the sphere "r0", on which the analysis will be
     * performed.  This loop must be executed before "pnt_nlat" is modified, as
     * it happens with the Driscoll--Healy grids below. */
    const REAL r0 = pnt->r[0];
    for (size_t i = 1; i < pnt_nlat; i++)
    {
        if (!CHARM(misc_is_nearly_equal)(pnt->r[i], r0, CHARM(glob_threshold)))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "All spherical radii in \"pnt->r\" must be "
                           "equal.");
            return;
        }
    }


    /* Get the maximum degree for which the grid in "pnt" was created. */
    const unsigned long nmax_grd =
                  CHARM(crd_point_quad_get_nmax_from_nlat)(pnt_type,
                                                           pnt->nlat);


    if (CHARM(crd_point_isDHGrid)(pnt_type) && (local_0_start == 0))
        /* See "shs_point_grd.c" for details on this */
        pnt_nlat += 1;


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


    /* Get the index of the equator for quadrature grids */
    const size_t equator = CHARM(crd_point_quad_equator)(pnt->type, nmax_grd);
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    /* Some useful constants */
    /* --------------------------------------------------------------------- */
    const size_t pnt_nlon     = CHARM(crd_point_get_local_nlon)(pnt);
    const size_t pnt_nlon_fft = pnt_nlon / 2 + 1;
    const _Bool even          = !(pnt_nlat % 2);
    const size_t nlatdo       = (pnt_nlat + 1 - even) / 2;
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);
    REAL c = PREC(0.0);
#if HAVE_MPI
    const size_t BLOCK_A = CHARM(glob_get_sha_block_lat_multiplier)();
#else
#   define BLOCK_A SIMD_BLOCK_A
#endif
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    REAL *r                      = NULL;
    REAL *ri                     = NULL;
    REAL *dm                     = NULL;
    INT *ips                     = NULL;
    REAL *ps                     = NULL;
    REAL *tv                     = NULL;
    REAL *uv                     = NULL;
    REAL *symmv                  = NULL;  /* "1.0" or "0.0" only */
    REAL *latsinv                = NULL;  /* "1.0" or "0.0" only */
    REAL *a                      = NULL;
    REAL *b                      = NULL;
    REAL *a2                     = NULL;
    REAL *b2                     = NULL;
    REAL *ftmp_in                = NULL;
    FFTWC(complex) *ftmp_out     = NULL;
    FFTW(plan) plan              = NULL;
    CHARM(shc_block) *shcs_block = NULL;
    MISC_SD_CALLOC_REAL_SIMD_INIT(t);
    MISC_SD_CALLOC_REAL_SIMD_INIT(u);
    MISC_SD_CALLOC_REAL_SIMD_INIT(symm);
    MISC_SD_CALLOC_REAL_SIMD_INIT(latsin);
    /* --------------------------------------------------------------------- */






    /* Initializations for recurrence relations to compute Legendre functions
     * */
    /* --------------------------------------------------------------------- */
    /* Prepare some variables to compute coefficients the "anm" and "bnm"
     * coefficients for the Legendre functions recurrences */
    r = (REAL *)calloc(2 * nmax + 4, sizeof(REAL));
    CHECK_NULL(r, BARRIER_1);


    ri = (REAL *)calloc(2 * nmax + 4, sizeof(REAL));
    CHECK_NULL(ri, BARRIER_1);


    CHARM(leg_func_r_ri)(nmax, r, ri);


    /* "dm" coefficients for sectorial Legendre functions */
    dm = (REAL *)calloc(nmax + 1, sizeof(REAL));
    CHECK_NULL(dm, BARRIER_1);


    CHARM(leg_func_dm)(nmax, r, ri, dm);
    /* --------------------------------------------------------------------- */






    /* Initialize block of spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    shcs_block = CHARM(shc_block_init)(shcs);
    CHECK_NULL(shcs_block, BARRIER_1);
    /* --------------------------------------------------------------------- */






    /* Auxiliary variable entering the computation of the lumped coefficients
     * */
    /* --------------------------------------------------------------------- */
    /* Here, we have to use the maximum degree "nmax_grd" for which the input
     * grid in "pnt" was created. */
    if ((pnt_type == CHARM_CRD_POINT_GRID_GL) ||
        (pnt_type == CHARM_CRD_POINT_GRID_DH1))
        c = PI / (REAL)(nmax_grd + 1);
    else if (pnt_type == CHARM_CRD_POINT_GRID_DH2)
        c = PI / (REAL)(2 * nmax_grd + 2);


    /* "4pi" normalization and normalization to "r0" and "shcs->mu" scaling
     * constants */
    c *= PREC(1.0) / (PREC(4.0) * PI) * (r0 / shcs->mu);
    /* --------------------------------------------------------------------- */






    /* Create a plan for FFT */
    /* --------------------------------------------------------------------- */
    {
#if HAVE_OPENMP && FFTW3_OMP
        if (FFTW(init_threads)() == 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFFTWINIT,
                           CHARM_ERR_FFTW_INIT_FAILURE);
            goto BARRIER_1;
        }


        FFTW(plan_with_nthreads)(omp_get_max_threads());
#endif


        REAL *x1           = NULL;
        FFTWC(complex) *x2 = NULL;


        x1 = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
        CHECK_NULL(x1, BARRIER_FFTW);


        x2 = (FFTWC(complex) *)FFTW(malloc)(pnt_nlon_fft *
                                            sizeof(FFTWC(complex)));
        CHECK_NULL(x2, BARRIER_FFTW);


        plan = FFTW(plan_dft_r2c_1d)(pnt_nlon, x1, x2, FFTW_ESTIMATE);
        if (plan == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);
            goto BARRIER_FFTW;
        }


BARRIER_FFTW:
        FFTW(free)(x1);
        FFTW(free)(x2);
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    /* Set all coefficients in "shcs" to zero */
    CHARM(shc_reset_coeffs)(shcs);


    MISC_SD_CALLOC_REAL_SIMD_ERR(t, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_ERR(u, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_ERR(symm, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_ERR(latsin, BLOCK_A, SIMD_BLOCK_A, err,
                                 BARRIER_1);
    const REAL_SIMD ROOT3_r = SET1_R(ROOT3);
#ifdef SIMD
    const RI_SIMD    zero_ri = SET_ZERO_RI;
    const RI_SIMD    one_ri  = SET1_RI(1);
    const RI_SIMD    mone_ri = SET1_RI(-1);
    const REAL_SIMD  zero_r  = SET_ZERO_R;
    const REAL_SIMD  BIG_r   = SET1_R(BIG);
    const REAL_SIMD  BIGI_r  = SET1_R(BIGI);
    const REAL_SIMD  BIGS_r  = SET1_R(BIGS);
    const REAL_SIMD  BIGSI_r = SET1_R(BIGSI);
    REAL_SIMD  tmp1_r,  tmp2_r;
    MASK_SIMD  mask1, mask2;
    MASK2_SIMD mask3;
#endif


    /* ................................................................. */
    ips = (INT *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                       SIMD_SIZE * BLOCK_A * nmax,
                                       sizeof(INT));
    CHECK_NULL(ips, BARRIER_1);


    ps = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                       SIMD_SIZE * BLOCK_A * nmax,
                                       sizeof(REAL));
    CHECK_NULL(ps, BARRIER_1);


    tv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE, sizeof(REAL));
    CHECK_NULL(tv, BARRIER_1);


    uv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE, sizeof(REAL));
    CHECK_NULL(uv, BARRIER_1);


    symmv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                          sizeof(REAL));
    CHECK_NULL(symmv, BARRIER_1);


    latsinv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
    CHECK_NULL(latsinv, BARRIER_1);


    a = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                      SIMD_SIZE * BLOCK_A *
                                      pnt_nlon_fft, sizeof(REAL));
    CHECK_NULL(a, BARRIER_1);


    b = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                      SIMD_SIZE * BLOCK_A *
                                      pnt_nlon_fft, sizeof(REAL));
    CHECK_NULL(b, BARRIER_1);


    a2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                       SIMD_SIZE * BLOCK_A *
                                       pnt_nlon_fft, sizeof(REAL));
    CHECK_NULL(a2, BARRIER_1);


    b2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                       SIMD_SIZE * BLOCK_A *
                                       pnt_nlon_fft, sizeof(REAL));
    CHECK_NULL(b2, BARRIER_1);


    ftmp_in = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
    CHECK_NULL(ftmp_in, BARRIER_1);


    ftmp_out = (FFTWC(complex) *)FFTW(malloc)(pnt_nlon_fft *
                                              sizeof(FFTWC(complex)));
    CHECK_NULL(ftmp_out, BARRIER_1);


BARRIER_1:
    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto FAILURE_1;
    /* ................................................................. */



    REAL_SIMD anms, bnms;
    REAL cw;
    REAL_SIMD wlf;
    _Bool npm_even; /* True if "n + m" is even */
    size_t l, ipv; /* "i + v" */
    int err_glob = 0;
    REAL_SIMD cs_sum;


    const size_t imax  = CHARM(shs_get_imax)(nlatdo, BLOCK_A, pnt);
    const size_t istep = SIMD_SIZE * BLOCK_A;


    /* Loop over latitudes */
    for (size_t i = 0; i < imax; i += istep)
    {
        for (l = 0; l < BLOCK_A; l++)
        {
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                /* Check whether the symmetry property of LFs needs to be
                 * applied */
                /* ----------------------------------------------------- */
                ipv = i + l * SIMD_SIZE + v;
                CHARM(crd_grd_check_symm)(ipv, v, local_0_start, equator,
                                          pnt_type, nlatdo, 0, even, symmv,
                                          latsinv);
                if (latsinv[v] == 1)
                {
                    tv[v] = SIN(pnt->lat[ipv]);
                    uv[v] = COS(pnt->lat[ipv]);
                }
                else
                {
                    tv[v] = uv[v] = PREC(0.0);
                    continue;
                }
                /* ----------------------------------------------------- */


                /* Lumped coefficients for the southern hemisphere
                 * (including the equator) */
                /* ----------------------------------------------------- */
                memcpy(ftmp_in, f + ipv * pnt_nlon, pnt_nlon * sizeof(REAL));
                FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


                cw = c * pnt->w[ipv];
                for (size_t j = 0; j < pnt_nlon_fft; j++)
                {
                    a[j * (SIMD_SIZE * BLOCK_A) +
                      l * SIMD_SIZE + v]  =  cw * ftmp_out[j][0];
                    b[j * (SIMD_SIZE * BLOCK_A) +
                      l * SIMD_SIZE + v]  = -cw * ftmp_out[j][1];
                }
                /* ----------------------------------------------------- */


                /* Lumped coefficients for the northern hemisphere */
                /* ----------------------------------------------------- */
                if (symmv[v])
                {
                    memcpy(ftmp_in, f + (pnt_nlat - ipv - 1) * pnt_nlon,
                           pnt_nlon * sizeof(REAL));
                    FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


                    cw = c * pnt->w[pnt_nlat - ipv - 1];
                    for (size_t j = 0; j < pnt_nlon_fft; j++)
                    {
                        a2[j * (SIMD_SIZE * BLOCK_A) +
                           l * SIMD_SIZE + v]  =  cw * ftmp_out[j][0];
                        b2[j * (SIMD_SIZE * BLOCK_A) +
                           l * SIMD_SIZE + v]  = -cw * ftmp_out[j][1];
                    }
                }
                /* ----------------------------------------------------- */
            }


            t[l]      = LOAD_R(&tv[0]);
            u[l]      = LOAD_R(&uv[0]);
            symm[l]   = LOAD_R(&symmv[0]);
            latsin[l] = LOAD_R(&latsinv[0]);


            /* Prepare arrays for sectorial Legendre functions */
            /* --------------------------------------------------------- */
            CHARM(leg_func_prepare)(uv, ps + l * SIMD_SIZE * nmax,
                                    ips + l * SIMD_SIZE * nmax, dm, nmax);
            /* --------------------------------------------------------- */
        }


        /* ------------------------------------------------------------- */
#if HAVE_OPENMP


#   undef SIMD_VARS1
#   ifdef SIMD
#       define SIMD_VARS1 shared(zero_ri, one_ri, mone_ri, zero_r) \
                     shared(BIG_r, BIGI_r, BIGS_r, BIGSI_r) \
                     private(tmp1_r, tmp2_r, mask1, mask2, mask3)
#   else
#       define SIMD_VARS1
#   endif


#   undef SIMD_VARS2
#   if HAVE_MPI
#       define SIMD_VARS2 shared(BLOCK_A)
#   else
#       define SIMD_VARS2
#   endif


#   define SIMD_VARS SIMD_VARS1 SIMD_VARS2


#pragma omp parallel default(none) \
shared(nmax, symm, r, ri, a, b, a2, b2, shcs, shcs_block, t, u, ps, ips) \
shared(latsin, pt, ROOT3_r, err_glob, err) \
private(anms, bnms, wlf) \
private(npm_even, l, cs_sum) SIMD_VARS
#endif
        {
        int err_priv = 0;


        REAL *anm = NULL;
        REAL *bnm = NULL;
        MISC_SD_CALLOC_REAL_SIMD_INIT(pnm0);
        MISC_SD_CALLOC_REAL_SIMD_INIT(pnm1);
        MISC_SD_CALLOC_REAL_SIMD_INIT(pnm2);
        MISC_SD_CALLOC_REAL_SIMD_INIT(x);
        MISC_SD_CALLOC_REAL_SIMD_INIT(y);
        MISC_SD_CALLOC_REAL_SIMD_INIT(z);
        MISC_SD_CALLOC_REAL_SIMD_INIT(amp);
        MISC_SD_CALLOC_REAL_SIMD_INIT(amm);
        MISC_SD_CALLOC_REAL_SIMD_INIT(bmp);
        MISC_SD_CALLOC_REAL_SIMD_INIT(bmm);
        MISC_SD_CALLOC_RI_SIMD_INIT(ix);
        MISC_SD_CALLOC_RI_SIMD_INIT(iy);
        MISC_SD_CALLOC_RI_SIMD_INIT(iz);
        MISC_SD_CALLOC_RI_SIMD_INIT(ixy);
        MISC_SD_CALLOC__BOOL_INIT(ds);


        MISC_SD_CALLOC_REAL_SIMD_ERR(pnm0, BLOCK_A, SIMD_BLOCK_A, err,
                                     BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(pnm1, BLOCK_A, SIMD_BLOCK_A, err,
                                     BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(pnm2, BLOCK_A, SIMD_BLOCK_A, err,
                                     BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(x, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(y, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(z, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(amp, BLOCK_A, SIMD_BLOCK_A, err,
                                     BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(amm, BLOCK_A, SIMD_BLOCK_A, err,
                                     BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(bmp, BLOCK_A, SIMD_BLOCK_A, err,
                                     BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(bmm, BLOCK_A, SIMD_BLOCK_A, err,
                                     BARRIER_2);
        MISC_SD_CALLOC_RI_SIMD_ERR(ix, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_2);
        MISC_SD_CALLOC_RI_SIMD_ERR(iy, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_2);
        MISC_SD_CALLOC_RI_SIMD_ERR(iz, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_2);
        MISC_SD_CALLOC_RI_SIMD_ERR(ixy, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_2);
        MISC_SD_CALLOC__BOOL_ERR(ds, BLOCK_A, SIMD_BLOCK_A, err, BARRIER_2);


        /* ............................................................. */
#if HAVE_MPI
        if (shcs_block->distributed)
        {
            CHARM(shc_block_reset_coeffs)(shcs_block);
            CHARM(shc_block_set_mfirst)(shcs_block, shcs, 0, err);
            if (!CHARM(err_isempty)(err))
            {
                err_priv = 1;
                goto BARRIER_2;
            }
        }
#endif
        /* ............................................................. */


        /* ............................................................. */
        anm = (REAL *)calloc(nmax + 1, sizeof(REAL));
        CHECK_NULL_OMP(anm, err_priv, BARRIER_2);


        bnm = (REAL *)calloc(nmax + 1, sizeof(REAL));
        CHECK_NULL_OMP(bnm, err_priv, BARRIER_2);


BARRIER_2:
        if (CHARM(err_omp_mpi)(&err_glob, &err_priv, CHARM_ERR_MALLOC_FAILURE,
                               CHARM_EMEM, err))
        {
#if HAVE_OPENMP
#pragma omp master
#endif
            if (!CHARM(err_isempty)(err))
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


#if HAVE_OPENMP
#pragma omp barrier
#endif
            goto FAILURE_2;
        }
        /* ............................................................. */


        /* Loop over harmonic orders */
        /* ------------------------------------------------------------- */
        /* For a more detailed description of the rational behind this block,
         * see "shs_point_sctr.c" */
        unsigned long m = shcs_block->mfirst;


        REAL_SIMD am, bm, a2m, b2m;


        do
        {
            /* The minimum and the maximum orders of the loop are the same for
             * all OpenMP threads */
            unsigned long mmin = shcs_block->mfirst;
            unsigned long mmax = CHARM_MIN(shcs_block->mlast, nmax);


#if HAVE_OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (m = mmin; m <= mmax; m++)
            {
                /* Apply polar optimization if asked to do so */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, &u[0],
                                                         BLOCK_A, pt))
                    continue;


                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Some useful substitutions */
                /* --------------------------------------------------------- */
                for (l = 0; l < BLOCK_A; l++)
                {
                    am  = LOAD_R(&a[(SIMD_SIZE * BLOCK_A) * m +
                                    l * SIMD_SIZE]);
                    bm  = LOAD_R(&b[(SIMD_SIZE * BLOCK_A) * m +
                                    l * SIMD_SIZE]);
                    a2m = LOAD_R(&a2[(SIMD_SIZE * BLOCK_A) * m +
                                     l * SIMD_SIZE]);
                    b2m = LOAD_R(&b2[(SIMD_SIZE * BLOCK_A) * m +
                                     l * SIMD_SIZE]);


                    amp[l] = MUL_R(ADD_R(am, MUL_R(symm[l], a2m)),
                                   latsin[l]);
                    amm[l] = MUL_R(SUB_R(am, MUL_R(symm[l], a2m)),
                                   latsin[l]);
                    bmp[l] = MUL_R(ADD_R(bm, MUL_R(symm[l], b2m)),
                                   latsin[l]);
                    bmm[l] = MUL_R(SUB_R(bm, MUL_R(symm[l], b2m)),
                                   latsin[l]);
                }
                /* --------------------------------------------------------- */


                unsigned long idx = CHARM(shc_block_get_idx)(shcs_block, m);


                /* Computation of spherical harmonic coefficients */
                if (m == 0)
                {
                    /* Zonal harmonics */
                    /* ----------------------------------------------------- */

                    /* P00 */
                    for (l = 0; l < BLOCK_A; l++)
                        pnm0[l] = SET1_R(PREC(1.0));
                    /* C00 */
                    CS_SUM(shcs_block->c[idx++], pnm0, amp);


                    if (nmax >= 1)
                    {
                        /* P10 */
                        for (l = 0; l < BLOCK_A; l++)
                            pnm1[l] = MUL_R(ROOT3_r, t[l]);
                        /* C10 */
                        CS_SUM(shcs_block->c[idx++], pnm1, amm);
                    }


                    /* P20, P30, ..., Pnmax,0 */
                    if (nmax >= 2)
                    {
                        /* Is "n + m" even?  Since we start the loop with "n
                         * = 2" and "m = 0", then the parity of the first "n
                         * + m" is always even.  Then, it changes with every
                         * loop iteration. */
                        npm_even = 1;


                        for (unsigned long n = 2; n <= nmax;
                             n++, npm_even = !npm_even)
                        {
                            anms = SET1_R(anm[n]);
                            bnms = SET1_R(bnm[n]);


                            for (l = 0; l < BLOCK_A; l++)
                            {
                                pnm2[l] = SUB_R(MUL_R(MUL_R(anms, t[l]),
                                                      pnm1[l]),
                                                MUL_R(bnms, pnm0[l]));
                                pnm0[l] = pnm1[l];
                                pnm1[l] = pnm2[l];
                            }


                            /* C20, C30, ..., Cnmax,0 */
                            CS_SUM(shcs_block->c[idx++], pnm2,
                                   npm_even ? amp : amm);
                        }
                    }
                    /* ----------------------------------------------------- */

                }
                else /* Non-zonal harmonics */
                {

                    /* Sectorial harmonics */
                    /* ----------------------------------------------------- */
                    for (l = 0; l < BLOCK_A; l++)
                    {
#ifdef SIMD
                        PNM_SECTORIAL_XNUM_SIMD(x[l], ix[l],
                                                ps[(SIMD_SIZE * nmax) * l +
                                                   (m - 1) * SIMD_SIZE],
                                                ips[(SIMD_SIZE * nmax) * l +
                                                   (m - 1) * SIMD_SIZE],
                                                pnm0[l],
                                                BIG_r, zero_r, zero_ri,
                                                mone_ri, mask1, mask2,
                                                SECTORIALS);
#else
                        PNM_SECTORIAL_XNUM(x[l], ix[l],
                                           ps[(SIMD_SIZE * nmax) * l +
                                              (m - 1) * SIMD_SIZE],
                                           ips[(SIMD_SIZE * nmax) * l +
                                               (m - 1) * SIMD_SIZE],
                                           pnm0[l]);
#endif
                    }


                    /* Cm,m; Sm,m */
                    CS_SUM(shcs_block->c[idx], pnm0, amp);
                    CS_SUM(shcs_block->s[idx++], pnm0, bmp);
                    /* ----------------------------------------------------- */


                    /* Tesseral harmonics */
                    /* ----------------------------------------------------- */
                    if (m < nmax)
                    {
                        anms = SET1_R(anm[m + 1]);
                        bnms = SET1_R(bnm[m + 1]);


                        for (l = 0; l < BLOCK_A; l++)
                        {
#ifdef SIMD
                            PNM_SEMISECTORIAL_XNUM_SIMD(x[l], y[l],
                                                        ix[l], iy[l],
                                                        wlf, t[l], anms,
                                                        pnm1[l],
                                                        mask1, mask2, mask3,
                                                        zero_r, zero_ri,
                                                        mone_ri,
                                                        BIG_r, BIGS_r,  BIGI_r,
                                                        SEMISECTORIALS);
#else
                            PNM_SEMISECTORIAL_XNUM(x[l], y[l], ix[l], iy[l],
                                                   wlf, t[l], anms, pnm1[l]);
#endif
                        }


                        /* Cm+1,m; Sm+1,m */
                        CS_SUM(shcs_block->c[idx], pnm1, amm);
                        CS_SUM(shcs_block->s[idx++], pnm1, bmm);


                        /* Loop over degrees */
                        /* ------------------------------------------------- */
                        for (l = 0; l < BLOCK_A; l++)
                            ds[l] = 0;


                        /* Is "n + m" even?  Since we start the loop with "n
                         * = m + 2", then the parity of the first "m + 2 + m"
                         * is always even.  Then, it changes with every loop
                         * iteration. */
                        npm_even = 1;


                        unsigned long n;
                        for (n = (m + 2);
                             CHARM(leg_func_use_xnum(ds, BLOCK_A)) &&
                             n <= nmax;
                             n++, npm_even = !npm_even)
                        {
                            anms = SET1_R(anm[n]);
                            bnms = SET1_R(bnm[n]);


                            /* Compute tesseral Legendre function */
                            for (l = 0; l < BLOCK_A; l++)
                            {
#ifdef SIMD
                                PNM_TESSERAL_XNUM_SIMD(x[l], y[l], z[l],
                                                       ix[l], iy[l], iz[l],
                                                       ixy[l], wlf, t[l],
                                                       anms, bnms,
                                                       pnm2[l], tmp1_r, tmp2_r,
                                                       mask1, mask2,
                                                       mask3, zero_r,
                                                       zero_ri, one_ri,
                                                       BIG_r, BIGI_r,
                                                       BIGS_r, BIGSI_r,
                                                       TESSERALS1, TESSERALS2,
                                                       ds[l]);
#else
                                PNM_TESSERAL_XNUM(x[l], y[l], z[l],
                                                  ix[l], iy[l], iz[l],
                                                  ixy[l], wlf, t[l],
                                                  anms, bnms, pnm2[l],
                                                  ds[l]);
#endif
                            }


                            /* Cm+2,m, Cm+3,m, ... and Sm+2,m, Sm+3,m, ... */
                            CS_SUM(shcs_block->c[idx], pnm2,
                                   npm_even ? amp : amm);
                            CS_SUM(shcs_block->s[idx++], pnm2,
                                   npm_even ? bmp : bmm);
                        }


                        if (n > nmax)
                            continue;


                        /* From now on, "F"-numbers can be used instead of the
                         * "X"-numbers to gain some speed */


                        /* We want to unroll the loop that follows.  To do
                         * that, we need to make sure that the loop starts with
                         * an even value of "n + m".  So if "n + m" is odd, we
                         * need to hard code one iteration. */
                        if (!npm_even)
                        {
                            LOOP_ITER(n, amm, bmm);
                            n++;
                        }


                        /* Now safely compute the rest of the loop (if any)
                         * using the "F"-numbers */
                        for (; (n + 1) <= nmax; n += 2)
                        {
                            LOOP_ITER(n,     amp, bmp);
                            LOOP_ITER(n + 1, amm, bmm);
                        }
                        /* ------------------------------------------------- */


                        if (n > nmax)
                            continue;


                        LOOP_ITER(n, amp, bmp);


                    } /* End of computation of tesseral harmonics */
                    /* ----------------------------------------------------- */


                } /* End of computation of spherical harmonic coefficients */
                /* --------------------------------------------------------- */


            } /* End of the loop over harmonic orders */


            /* See "shs_point_sctr.c" for the rational here */
            m = mmax + 1;


#if HAVE_MPI
            /* Now that we have processed all harmonic degree "mmin", ...,
             * "mmax", it is time to sum the contributions from the MPI
             * processes and send the result to the processes holding this
             * chunk of coefficients */
            CHARM(shc_block_set_coeffs)(shcs, shcs_block, mmin, mmax, err);
            /* Only the master thread in "shc_block_set_coeffs" modifies "err",
             * so we do not need below "CHARM(err_omp_mpi)" */
            if (CHARM(err_omp_mpi)(&err_glob, &err_priv,
                                   CHARM_ERR_MALLOC_FAILURE, CHARM_EMEM, err))
            {
#   if HAVE_OPENMP
#pragma omp master
#   endif
                if (!CHARM(err_isempty)(err))
                    CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


#   if HAVE_OPENMP
#pragma omp barrier
#   endif
                goto FAILURE_2;
            }
#endif
        }
        while (m <= nmax);
        /* ------------------------------------------------------------- */


FAILURE_2:
        free(anm);
        free(bnm);
        MISC_SD_FREE(pnm0);
        MISC_SD_FREE(pnm1);
        MISC_SD_FREE(pnm2);
        MISC_SD_FREE(x);
        MISC_SD_FREE(y);
        MISC_SD_FREE(z);
        MISC_SD_FREE(amp);
        MISC_SD_FREE(amm);
        MISC_SD_FREE(bmp);
        MISC_SD_FREE(bmm);
        MISC_SD_FREE(ix);
        MISC_SD_FREE(iy);
        MISC_SD_FREE(iz);
        MISC_SD_FREE(ixy);
        MISC_SD_FREE(ds);
        }
        /* ------------------------------------------------------------- */


    } /* End of the loop over latitude parallels */
    /* ----------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE_1:
    free(r);
    free(ri);
    free(dm);
    MISC_SD_FREE(t);
    MISC_SD_FREE(u);
    MISC_SD_FREE(symm);
    MISC_SD_FREE(latsin);
    CHARM(free_aligned)(ips);
    CHARM(free_aligned)(ps);
    CHARM(free_aligned)(tv);
    CHARM(free_aligned)(uv);
    CHARM(free_aligned)(symmv);
    CHARM(free_aligned)(latsinv);
    CHARM(free_aligned)(a);
    CHARM(free_aligned)(b);
    CHARM(free_aligned)(a2);
    CHARM(free_aligned)(b2);
    FFTW(free)(ftmp_in);
    FFTW(free)(ftmp_out);
    FFTW(destroy_plan)(plan);
#if HAVE_OPENMP && FFTW3_OMP
    FFTW(cleanup_threads)();
#else
    FFTW(cleanup)();
#endif
    CHARM(shc_block_free)(shcs_block);


    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto EXIT;
    /* --------------------------------------------------------------------- */






    /* Rescale the coefficients to the sphere of radius "shcs->r0" */
    /* --------------------------------------------------------------------- */
    if (!CHARM(misc_is_nearly_equal)(shcs->r, r0, CHARM(glob_threshold)))
    {
        /* This is the desired radius to scale the output coefficients */
        const REAL rtmp = shcs->r;


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
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
EXIT:
#if HAVE_MPI
    CHARM(mpi_err_gather)(err);
#endif


    return;
    /* --------------------------------------------------------------------- */
}

