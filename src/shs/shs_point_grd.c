/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#if HAVE_OPENMP
#   include <omp.h>
#endif
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"
#include "../crd/crd_point_isDHGrid.h"
#include "../crd/crd_point_quad_equator.h"
#include "../crd/crd_point_quad_get_nmax_from_nlat.h"
#include "../shc/shc_block_struct.h"
#include "../shc/shc_block_init.h"
#include "../shc/shc_block_free.h"
#include "../shc/shc_block_have_order.h"
#include "../shc/shc_block_get_coeffs.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_enm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_isempty_all_mpi_processes.h"
#include "../err/err_omp_mpi.h"
#include "../crd/crd_grd_check_symm.h"
#include "../crd/crd_point_isQuadGrid.h"
#include "../crd/crd_point_isCustGrid.h"
#include "../crd/crd_point_get_local_nlat.h"
#include "../crd/crd_point_get_local_nlon.h"
#include "../crd/crd_point_get_local_0_start.h"
#include "../misc/misc_arr_chck_lin_incr.h"
#include "../misc/misc_arr_chck_symm.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#include "../misc/misc_sd_calloc.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
#include "../glob/glob_get_shs_block_lat_multiplier.h"
#include "shs_grd_fft.h"
#include "shs_grd_point_fft_check.h"
#include "shs_grd_fft_lc.h"
#include "shs_grd_lr.h"
#include "shs_grd_lr2.h"
#include "shs_point_kernels.h"
#include "shs_r_eq_rref.h"
#include "shs_get_mur_dorder_npar.h"
#include "shs_point_gradn.h"
#include "shs_lc_struct.h"
#include "shs_lc_init.h"
#include "shs_lc_free.h"
#include "shs_max_npar.h"
#include "shs_ratios.h"
#include "shs_get_imax.h"
#include "shs_point_grd.h"
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
/* ------------------------------------------------------------------------- */






void CHARM(shs_point_grd)(const CHARM(point) *pnt,
                          const CHARM(shc) *shcs,
                          unsigned long nmax,
                          int dr,
                          int dlat,
                          int dlon,
                          REAL **f,
                          CHARM(err) *err)
{
    /* Check the latitudes */
    /* --------------------------------------------------------------------- */

    /* Check whether the number of latitudes is even or odd */
    /* ..................................................................... */
    const int pnt_type = pnt->type;
    size_t pnt_nlat = CHARM(crd_point_get_local_nlat)(pnt);
    const size_t local_0_start = CHARM(crd_point_get_local_0_start)(pnt);


    /* For the full Driscoll--Healy grids, "pnt_nlat" is always an even number
     * and the grid is not symmetric in terms of our definition (the north pole
     * does not have its negative counterpart -- the south pole).  However, we
     * know that except for the north pole, the Driscoll--Healy grids *are*
     * symmetric, so the symmetry property of Legendre functions could be used
     * if the north pole is treated properly.  To this end, let's increase the
     * number of points in "pnt_nlat", so that we can apply our algorithm for
     * symmetric grids.  After increasing "pnt_nlat", its value is odd, so set
     * "even" to zero.
     *
     * With MPI, we have to increase "pnt_nlat" only for the chunk that has
     * the north pole, which is always the one with "local_0_start == 0" due to
     * our restrictions. */
    if (CHARM(crd_point_isDHGrid)(pnt_type) && (local_0_start == 0))
        pnt_nlat += 1;


    const _Bool even = !(pnt_nlat % 2);
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
    else if (CHARM(crd_point_isQuadGrid)(pnt_type))
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


    /* Get the index of the equator for quadrature grids */
    /* ..................................................................... */
    const unsigned long nmax_grd =
                CHARM(crd_point_quad_get_nmax_from_nlat)(pnt->type, pnt->nlat);
    const size_t equator = CHARM(crd_point_quad_equator)(pnt->type, nmax_grd);
    /* ..................................................................... */


    /* Finally, if the grid is symmetric, we modify the number of latitudes
     * "pnt_nlat" to be equal to the number of latitudes on one hemisphere
     * only (including the equator if present). This is because the
     * time-consuming "for loop" over evaluation points now needs to run for
     * one hemisphere only, while the results for the other hemisphere are
     * obtained by exploiting the symmetry property of Legendre functions. This
     * reduces the number of Legendre functions that need to be evaluated by
     * a factor of ~2, so saves some computational time */
    const size_t nlatdo = (symm) ? (pnt_nlat + 1 - even) / 2 : pnt_nlat;
    /* --------------------------------------------------------------------- */






    /* Check the longitudes.  If possible, FFT is employed along the
     * latitudinal parallels.  Otherwise, the PSLR algorithm is used.  The
     * latter is slower, but can be used for any grid with a constant
     * longitudinal sampling.  Below, we determined whether FFT can be applied
     * or not. */
    /* --------------------------------------------------------------------- */
    const size_t pnt_nlon = CHARM(crd_point_get_local_nlon)(pnt);


    /* If "pnt" is a user-defined grid with more than one longitude, we have to
     * check whether the longitudinal step is constant.  For quadrature grids,
     * the longitudinal step is constant by definition, so no check is
     * needed. */
    if (CHARM(crd_point_isCustGrid)(pnt->type) && (pnt_nlon > 1))
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
    const REAL deltalon = (pnt_nlon > 1) ? pnt->lon[1] - pnt->lon[0] :
                                           PREC(0.0);


    /* Auxiliary constant to be used only in case the PSLR algorithm is
     * applied.  To suppress, a compiler warning, the value is initialized to
     * zero. */
    REAL lon0 = PREC(0.0);


    /* Length of the lumped coefficients arrays in case FFT will be applied */
    const size_t nfc = pnt_nlon / 2 + 1;


    const _Bool use_fft = CHARM(shs_grd_point_fft_check)(pnt, deltalon, nmax);
    if (!use_fft)
    {
        /* Get the origin of the longitude "pnt->lon" vector (will be necessary
         * later for the PSLR algorithm) */
        lon0 = pnt->lon[0];
    }
    /* --------------------------------------------------------------------- */






    /* Check whether all values of "pnt->r" are equal to "shcs->r".  If true,
     * a faster code can be used inside "shs_point_kernel".  */
    /* --------------------------------------------------------------------- */
    const _Bool r_eq_rref = CHARM(shs_r_eq_rref)(pnt, shcs);
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    REAL *r                      = NULL;
    REAL *ri                     = NULL;
    REAL *dm                     = NULL;
    FFTW(plan) plan              = NULL;
    INT *ips                     = NULL;
    REAL *ps                     = NULL;
    REAL *latv                   = NULL;
    REAL *tv                     = NULL;
    REAL *uv                     = NULL;
    REAL *symmv                  = NULL;
    REAL *latsinv                = NULL;
    REAL *pnt_rv                 = NULL;
    REAL *pnt_r2v                = NULL;
    FFTWC(complex) *fc           = NULL;
    FFTWC(complex) *fc2          = NULL;
    REAL *ftmp                   = NULL;
    REAL *fi                     = NULL;
    REAL *fi2                    = NULL;
    REAL *fc_simd                = NULL;
    REAL *fc2_simd               = NULL;
    CHARM(shc_block) *shcs_block = NULL;
    MISC_SD_CALLOC_REAL_SIMD_INIT(t);
    MISC_SD_CALLOC_REAL_SIMD_INIT(u);
    MISC_SD_CALLOC_REAL_SIMD_INIT(symm_simd);
    MISC_SD_CALLOC_REAL_SIMD_INIT(ratio);
    MISC_SD_CALLOC_REAL_SIMD_INIT(ratio2);
    /* --------------------------------------------------------------------- */






    /* Get some constants */
    /* --------------------------------------------------------------------- */
    REAL mur;  /* "(shcs->mu / shcs->r)^dorder" */
    unsigned dorder;  /* "0" for potential,
                         "1" for first-order derivatives,
                         "2" for second-order derivatives */
    size_t npar; /* Number of quantities to be synthesized:
                    "3" for first-order derivatives,
                    "6" for second-order derivatives,
                    "1" otherwise. */
    CHARM(shs_get_mur_dorder_npar)(shcs, dr, dlat, dlon, &mur, &dorder, &npar,
                                   err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto FAILURE_1;


    int grad;
    if ((dr == GRAD_1) && (dlat == GRAD_1) && (dlon == GRAD_1))
        grad = 1;
    else if ((dr == GRAD_2) && (dlat == GRAD_2) && (dlon == GRAD_2))
        grad = 2;
    else
        grad = 0;


#if HAVE_MPI
    const size_t BLOCK_S = CHARM(glob_get_shs_block_lat_multiplier)();
#else
#   define BLOCK_S SIMD_BLOCK_S
#endif
    const size_t nfi_1par = pnt_nlon * SIMD_SIZE * BLOCK_S;
    const size_t nfi = npar * nfi_1par;
    /* --------------------------------------------------------------------- */






    /* Initializations for recurrence relations to compute Legendre
     * functions. */
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






    /* Create a FFT plan */
    /* --------------------------------------------------------------------- */
    if (use_fft)
    {
#if HAVE_OPENMP && FFTW3_OMP
        if (FFTW(init_threads)() == 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFFTWINIT,
                           "FFTW failed to initialize threads.");
            goto BARRIER_1;
        }


        FFTW(plan_with_nthreads)(omp_get_max_threads());
#endif


        REAL *x1           = NULL;
        FFTWC(complex) *x2 = NULL;


        x1 = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
        CHECK_NULL(x1, BARRIER_FFTW);


        x2 = (FFTWC(complex) *)FFTW(malloc)(nfc * sizeof(FFTWC(complex)));
        CHECK_NULL(x2, BARRIER_FFTW);


        plan = FFTW(plan_dft_c2r_1d)(pnt_nlon, x2, x1, FFTW_ESTIMATE);
        if (plan == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);
            goto BARRIER_FFTW;
        }


BARRIER_FFTW:
        FFTW(free)(x1);
        FFTW(free)(x2);


        if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
            goto BARRIER_1;
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    ips = (INT *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                       nmax * SIMD_SIZE * BLOCK_S,
                                       sizeof(INT));
    CHECK_NULL(ips, BARRIER_1);


    ps = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                       nmax * SIMD_SIZE * BLOCK_S,
                                       sizeof(REAL));
    CHECK_NULL(ps, BARRIER_1);


    latv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                         sizeof(REAL));
    CHECK_NULL(latv, BARRIER_1);


    tv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE, sizeof(REAL));
    CHECK_NULL(tv, BARRIER_1);


    uv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE, sizeof(REAL));
    CHECK_NULL(uv, BARRIER_1);


    symmv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                          SIMD_SIZE * BLOCK_S,
                                          sizeof(REAL));
    CHECK_NULL(symmv, BARRIER_1);


    latsinv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                            SIMD_SIZE * BLOCK_S,
                                            sizeof(REAL));
    CHECK_NULL(latsinv, BARRIER_1);


    pnt_rv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                           sizeof(REAL));
    CHECK_NULL(pnt_rv, BARRIER_1);


    if (use_fft)
    {
        ftmp = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
        CHECK_NULL(ftmp, BARRIER_1);


        fc = (FFTWC(complex) *)FFTW(malloc)(nfc * sizeof(FFTWC(complex)));
        CHECK_NULL(fc, BARRIER_1);
        memset(fc, 0, nfc * sizeof(FFTWC(complex)));


        fc_simd = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                npar * nfc * SIMD_SIZE *
                                                BLOCK_S * 2,
                                                sizeof(REAL));
        CHECK_NULL(fc_simd, BARRIER_1);
    }
    else
    {
        fi = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nfi,
                                           sizeof(REAL));
        CHECK_NULL(fi, BARRIER_1);
    }


    if (symm)
    {
        pnt_r2v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                SIMD_SIZE, sizeof(REAL));
        CHECK_NULL(pnt_r2v, BARRIER_1);


        if (use_fft)
        {
            fc2 = (FFTWC(complex) *)FFTW(malloc)(nfc * sizeof(FFTWC(complex)));
            CHECK_NULL(fc2, BARRIER_1);
            memset(fc2, 0, nfc * sizeof(FFTWC(complex)));


            fc2_simd = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                     npar * nfc * SIMD_SIZE *
                                                     BLOCK_S * 2,
                                                     sizeof(REAL));
            CHECK_NULL(fc2_simd, BARRIER_1);
        }
        else
        {
            fi2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nfi,
                                                sizeof(REAL));
            CHECK_NULL(fi2, BARRIER_1);
        }
    }


    shcs_block = CHARM(shc_block_init)(shcs);
    CHECK_NULL(shcs_block, BARRIER_1);


BARRIER_1:
    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto FAILURE_1;
    /* --------------------------------------------------------------------- */






    /* Loop over grid latitudes */
    /* ----------------------------------------------------------------- */
    size_t ipv, l;
    int err_glob = 0;
    unsigned long lc_err_glob = 0;
    const size_t size_blk2 = SIMD_SIZE * BLOCK_S * 2;


    MISC_SD_CALLOC_REAL_SIMD_ERR(t, BLOCK_S, SIMD_BLOCK_S, err, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_ERR(u, BLOCK_S, SIMD_BLOCK_S, err, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_ERR(symm_simd, BLOCK_S, SIMD_BLOCK_S, err,
                                 BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_ERR(ratio, BLOCK_S, SIMD_BLOCK_S, err, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_ERR(ratio2, BLOCK_S, SIMD_BLOCK_S, err,
                                 BARRIER_1);
    REAL_SIMD pnt_r, pnt_r2;
    for (l = 0; l < BLOCK_S; l++)
        ratio[l] = ratio2[l] = SET_ZERO_R;


    /* Get the polar optimization threshold */
    const REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


    /* Radius of the reference sphere that is associated with the spherical
     * harmonic coefficients */
    const REAL_SIMD rref = SET1_R(shcs->r);


    const size_t imax  = CHARM(shs_get_imax)(nlatdo, BLOCK_S, pnt);
    const size_t istep = SIMD_SIZE * BLOCK_S;


    for (size_t i = 0; i < imax; i += istep)
    {
        for (l = 0; l < BLOCK_S; l++)
        {
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                /* Check whether the symmetry property of LFs can be
                 * applied */
                /* ----------------------------------------------------- */
                ipv = i + l * SIMD_SIZE + v;
                CHARM(crd_grd_check_symm)(ipv, v, local_0_start, equator,
                                          pnt_type,
                                          nlatdo, symm, even,
                                          symmv + l * SIMD_SIZE,
                                          latsinv + l * SIMD_SIZE);


                if (latsinv[l * SIMD_SIZE + v] == 1)
                {
                    latv[v]   = pnt->lat[ipv];
                    tv[v]     = SIN(pnt->lat[ipv]);
                    uv[v]     = COS(pnt->lat[ipv]);
                    pnt_rv[v] = pnt->r[ipv];


                    if (symm)
                    {
                        if (CHARM(crd_point_isDHGrid)(pnt_type) &&
                            ((ipv + local_0_start) == 0))
                            /* For the Driscoll--Healy grids, we increased
                             * "pnt_nlat" by one, so to ensure we won't
                             * read outside the bounds of "pnt->r", we set
                             * now the non-existing radius at the south
                             * pole to zero. */
                            pnt_r2v[v] = PREC(0.0);
                        else
                            pnt_r2v[v] = pnt->r[pnt_nlat - ipv - 1];
                    }
                }
                else
                {
                    latv[v] = tv[v] = uv[v] = pnt_rv[v] = PREC(0.0);
                    if (symm)
                        pnt_r2v[v] = PREC(0.0);


                    continue;
                }
                /* ----------------------------------------------------- */
            }


            t[l]         = LOAD_R(&tv[0]);
            u[l]         = LOAD_R(&uv[0]);
            pnt_r        = LOAD_R(&pnt_rv[0]);
            symm_simd[l] = LOAD_R(&symmv[l * SIMD_SIZE]);


            ratio[l]  = DIV_R(rref, pnt_r);
            if (symm)
            {
                pnt_r2    = LOAD_R(&pnt_r2v[0]);
                ratio2[l] = DIV_R(rref, pnt_r2);
            }


            /* Prepare arrays for sectorial Legendre functions */
            /* --------------------------------------------------------- */
            CHARM(leg_func_prepare)(uv, ps + l * SIMD_SIZE * nmax,
                                    ips + l * SIMD_SIZE * nmax, dm, nmax);
            /* --------------------------------------------------------- */


            if (use_fft)
            {
                /* Reset the lumped coefficients.  Required in some
                 * cases. */
                /* ----------------------------------------------------- */
                memset(fc, 0, nfc * sizeof(FFTWC(complex)));


                if (symm)
                    memset(fc2, 0, nfc * sizeof(FFTWC(complex)));
                /* ----------------------------------------------------- */
            }
            else
            {
                /* The "fi" vector represents the synthesized quantity "f"
                 * for the "ipv"th latitude parallel. Therefore, it needs
                 * to be reinitialized to zero for each "ith" latitude. The
                 * same holds true for "fi2" in case of symmetric grids. */
                /* ----------------------------------------------------- */
                memset(fi, 0, nfi * sizeof(REAL));


                if (symm)
                    memset(fi2, 0, nfi * sizeof(REAL));
                /* ----------------------------------------------------- */
            }
        }


        /* ------------------------------------------------------------- */
#undef MPI_VARS
#if HAVE_MPI
        _Bool have_order;
#   define MPI_VARS shared(have_order, BLOCK_S)
#else
#   define MPI_VARS
#endif


#if HAVE_OPENMP
#pragma omp parallel default(none) \
shared(nmax, err, dorder, pt, t, u, ri, r, ips, ps, dr, dlat, dlon) \
shared(symm_simd, ratio, ratio2, r_eq_rref, shcs, shcs_block) \
shared(use_fft, nfc, symm, grad, fi, fi2, nfi, nfi_1par) \
shared(pnt_type, pnt_nlon, deltalon, lon0, fc_simd, fc2_simd, err_glob) \
shared(lc_err_glob) \
private(l) MPI_VARS
#endif
        {
        /* ------------------------------------------------------------- */
        int err_priv = 0;
        unsigned long lc_err_priv = 0;


        CHARM(lc) *lc    = NULL;
        REAL *anm        = NULL;
        REAL *bnm        = NULL;
        REAL *enm        = NULL;
        REAL *fi_thread  = NULL;
        REAL *fi2_thread = NULL;
        MISC_SD_CALLOC_REAL_SIMD_INIT(ratiom);
        MISC_SD_CALLOC_REAL_SIMD_INIT(ratio2m);
        /* ------------------------------------------------------------- */


        /* Feed "shcs_block" starting with order "0" */
        /* ------------------------------------------------------------- */
        CHARM(shc_block_get_coeffs)(shcs
#if HAVE_MPI
                                    , shcs_block,
                                    0,
                                    err
#endif
                                   );
        if (!CHARM(err_isempty)(err))
        {
#if HAVE_OPENMP
#pragma omp master
#endif
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


#if HAVE_OPENMP
#pragma omp barrier
#endif
            goto BARRIER_2;
        }
        /* ------------------------------------------------------------- */


        /* ------------------------------------------------------------- */
        lc = CHARM(shs_lc_init)();
        CHECK_NULL_OMP(lc, err_priv, BARRIER_2);


        anm = (REAL *)calloc(nmax + 1, sizeof(REAL));
        CHECK_NULL_OMP(anm, err_priv, BARRIER_2);


        bnm = (REAL *)calloc(nmax + 1, sizeof(REAL));
        CHECK_NULL_OMP(bnm, err_priv, BARRIER_2);


        if (dorder > 0)
        {
            enm = (REAL *)calloc(nmax + 1, sizeof(REAL));
            CHECK_NULL_OMP(enm, err_priv, BARRIER_2);
        }



        if (!use_fft)
        {
            fi_thread = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nfi,
                                                      sizeof(REAL));
            CHECK_NULL_OMP(fi_thread, err_priv, BARRIER_2);


            fi2_thread = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nfi,
                                                       sizeof(REAL));
            CHECK_NULL_OMP(fi2_thread, err_priv, BARRIER_2);
        }


        MISC_SD_CALLOC_REAL_SIMD_ERR(ratiom, BLOCK_S, SIMD_BLOCK_S, err,
                                     BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(ratio2m, BLOCK_S, SIMD_BLOCK_S, err,
                                     BARRIER_2);
        for (l = 0; l < BLOCK_S; l++)
            ratiom[l] = ratio2m[l] = SET_ZERO_R;
        /* ------------------------------------------------------------- */


        /* Loop over harmonic orders */
        /* ------------------------------------------------------------- */
        /* For a more detailed description of the rational behind this block,
         * see "shs_point_sctr.c" */
        unsigned long m = shcs_block->mfirst;


        do
        {
#if HAVE_MPI
#if HAVE_OPENMP
#pragma omp master
#endif
            have_order = CHARM(shc_block_have_order)(shcs_block, m);


#if HAVE_OPENMP
#pragma omp barrier
#endif

            if (shcs->distributed && !have_order)
            {
                CHARM(shc_block_get_coeffs)(shcs, shcs_block, m, err);
                if (!CHARM(err_isempty)(err))
                {
                    CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                    goto BARRIER_2;
                }
            }
#endif


            /* ............................................................. */
BARRIER_2:
            if (CHARM(err_omp_mpi)(&err_glob, &err_priv,
                                   CHARM_ERR_MALLOC_FAILURE, CHARM_EMEM, err))
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
                                                         BLOCK_S, pt))
                    continue;


                /* Compute "(R / r)^(m + 1)" ("m" is not a typo) */
                CHARM(shs_ratios)(ratio, ratio2, symm, m, ratiom, ratio2m);


                /* "anm" and "bnm" coefficients for Legendre recurrence
                 * relations and for their derivatives */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);
                if (dorder > 0)
                    CHARM(leg_func_enm)(nmax, m, r, ri, enm);


                /* Computation of the lumped coefficients */
                /* --------------------------------------------------------- */
#undef KERNEL_IO_PARS
#define KERNEL_IO_PARS (nmax, m, shcs_block, r_eq_rref, anm, bnm, enm,        \
                        &t[0], &u[0], ps, ips,                                \
                        &ratio[0], &ratio2[0], &ratiom[0], &ratio2m[0],       \
                        &symm_simd[0], dorder, lc);
                if ((dr == 0) && (dlat == 0) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr0_dlat0_dlon0) KERNEL_IO_PARS;
                }
                else if ((dr == 1) && (dlat == 0) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr1_dlat0_dlon0) KERNEL_IO_PARS;
                }
                else if ((dr == 2) && (dlat == 0) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr2_dlat0_dlon0) KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 1) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr0_dlat1_dlon0) KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 2) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr0_dlat2_dlon0) KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 0) && (dlon == 1))
                {
                    CHARM(shs_point_kernel_dr0_dlat0_dlon1) KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 0) && (dlon == 2))
                {
                    CHARM(shs_point_kernel_dr0_dlat0_dlon2) KERNEL_IO_PARS;
                }
                else if ((dr == 1) && (dlat == 1) && (dlon == 0))
                {
                    CHARM(shs_point_kernel_dr1_dlat1_dlon0) KERNEL_IO_PARS;
                }
                else if ((dr == 1) && (dlat == 0) && (dlon == 1))
                {
                    CHARM(shs_point_kernel_dr1_dlat0_dlon1) KERNEL_IO_PARS;
                }
                else if ((dr == 0) && (dlat == 1) && (dlon == 1))
                {
                    CHARM(shs_point_kernel_dr0_dlat1_dlon1) KERNEL_IO_PARS;
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


                if (lc->error)
                    lc_err_priv += 1;
                /* --------------------------------------------------------- */


                /* --------------------------------------------------------- */
                if (use_fft)
                    CHARM(shs_grd_fft_lc)(m, deltalon, grad, lc,
                                          symm, &symm_simd[0], pnt_type,
                                          nfc, fc_simd, fc2_simd);
                else
                    CHARM(shs_grd_lr)(m, lon0, deltalon, pnt_nlon, pnt_type,
                                      grad, nfi_1par, lc, symm,
                                      fi_thread, fi2_thread);


                if (lc->error)
                    lc_err_priv += 1;
                /* --------------------------------------------------------- */


            } /* End of the loop over harmonic orders */


            /* See "shs_point_sctr.c" for the rational here */
            m = mmax + 1;
        }
        while (m <= nmax);


        /* With the Chebyshev recurrences, we must gather the results from
         * threads and combine them */
        if (!use_fft)
        {
#if HAVE_OPENMP
#pragma omp critical
#endif
            {
            for (l = 0; l < nfi; l++)
                fi[l] += fi_thread[l];


            lc_err_glob += lc_err_priv;
            }

            if (symm)
            {
#if HAVE_OPENMP
#pragma omp critical
#endif
                for (l = 0; l < nfi; l++)
                    fi2[l] += fi2_thread[l];
            }
        }


FAILURE_2:
        CHARM(shs_lc_free)(lc);
        free(anm);
        free(bnm);
        free(enm);
        CHARM(free_aligned)(fi_thread);
        CHARM(free_aligned)(fi2_thread);
        MISC_SD_FREE(ratiom);
        MISC_SD_FREE(ratio2m);
        }  /* End of parallel block */
        /* ------------------------------------------------------------- */


        for (size_t p = 0; p < npar; p++)
        {
            if (use_fft)
                /* Fourier transform along the latitude parallels */
                CHARM(shs_grd_fft)(i, pnt_type, pnt_nlat, pnt_nlon,
                                   latsinv, NULL, NULL, PREC(0.0),
                                   fc, fc2, nfc,
                                   &fc_simd[p * nfc * size_blk2],
                                   &fc2_simd[p * nfc * size_blk2],
                                   mur, plan, symmv,
                                   ftmp, f[p]);
            else
                CHARM(shs_grd_lr2)(i, latsinv, pnt_type, pnt_nlat,
                                   pnt_nlon, symmv, mur, NULL, NULL,
                                   PREC(0.0),
                                   &fi[p * nfi_1par], &fi2[p * nfi_1par],
                                   f[p]);
        }


    } /* End of the loop over latitude parallels */


    if (lc_err_glob)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
    /* ----------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE_1:
    free(r);
    free(ri);
    free(dm);
    if (use_fft)
    {
        FFTW(destroy_plan)(plan);
#if HAVE_OPENMP && FFTW3_OMP
        FFTW(cleanup_threads)();
#else
        FFTW(cleanup)();
#endif
    }
    CHARM(free_aligned)(ips);
    CHARM(free_aligned)(ps);
    CHARM(free_aligned)(latv);
    CHARM(free_aligned)(tv);
    CHARM(free_aligned)(uv);
    CHARM(free_aligned)(symmv);
    CHARM(free_aligned)(latsinv);
    CHARM(free_aligned)(pnt_rv);
    FFTW(free)(ftmp);
    FFTW(free)(fc);
    CHARM(free_aligned)(fi);
    CHARM(free_aligned)(pnt_r2v);
    FFTW(free)(fc2);
    CHARM(free_aligned)(fi2);
    CHARM(free_aligned)(fc_simd);
    CHARM(free_aligned)(fc2_simd);
    CHARM(shc_block_free)(shcs_block);
    MISC_SD_FREE(t);
    MISC_SD_FREE(u);
    MISC_SD_FREE(symm_simd);
    MISC_SD_FREE(ratio);
    MISC_SD_FREE(ratio2);
    /* --------------------------------------------------------------------- */






    return;
}
