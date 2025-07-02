/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"
#include "../shc/shc_block_struct.h"
#include "../shc/shc_block_init.h"
#include "../shc/shc_block_free.h"
#include "../shc/shc_block_have_order.h"
#include "../shc/shc_block_get_coeffs.h"
#include "../crd/crd_point_get_local_nlat.h"
#include "../crd/crd_point_get_local_nlon.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_enm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#include "../misc/misc_sd_calloc.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_isempty_all_mpi_processes.h"
#include "../err/err_omp_mpi.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
#include "../glob/glob_get_shs_block_lat_multiplier.h"
#include "shs_point_kernels.h"
#include "shs_sctr_mulc.h"
#include "shs_r_eq_rref.h"
#include "shs_lc_struct.h"
#include "shs_lc_init.h"
#include "shs_lc_free.h"
#include "shs_get_mur_dorder_npar.h"
#include "shs_point_gradn.h"
#include "shs_max_npar.h"
#include "shs_get_imax.h"
#include "shs_rpows.h"
#include "shs_point_sctr.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
#undef LC_CS
#define LC_CS(x)                                                              \
        fi_thread[idx] = ADD_R(fi_thread[idx],                                \
                               ADD_R(MUL_R(lc->CAT(a, x)[l], clonim),         \
                                     MUL_R(lc->CAT(b, x)[l], slonim)));


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
    char err_msg[CHARM_ERR_MAX_MSG];
    REAL *r                      = NULL;
    REAL *ri                     = NULL;
    REAL *dm                     = NULL;
    CHARM(shc_block) *shcs_block = NULL;
    INT *ips                     = NULL;
    REAL *ps                     = NULL;
    REAL *tv                     = NULL;
    REAL *uv                     = NULL;
    REAL *lonv                   = NULL;
    REAL *pnt_rv                 = NULL;
    REAL *tmpv                   = NULL;
    REAL_SIMD *rpows             = NULL;
    MISC_SD_CALLOC_REAL_SIMD_INIT(t);
    MISC_SD_CALLOC_REAL_SIMD_INIT(u);
    MISC_SD_CALLOC_REAL_SIMD_INIT(fi);


    int grad = -1;


    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


    /* Radius of the reference sphere that is associated with the spherical
     * harmonic coefficients */
    REAL_SIMD rref = SET1_R(shcs->r);


    /* Check whether all values of "pnt->r" are equal to "shcs->r".  If true,
     * a faster code can be used inside "shs_point_kernel".  */
    _Bool r_eq_rref = CHARM(shs_r_eq_rref)(pnt, shcs);


#if HAVE_MPI
    size_t BLOCK_S = CHARM(glob_get_shs_block_lat_multiplier)();
#else
#   define BLOCK_S SIMD_BLOCK_S
#endif
    /* --------------------------------------------------------------------- */


    /* The number of latitudes and longitudes must match */
    /* --------------------------------------------------------------------- */
    /* Get the number of latitudes stored locally */
    size_t local_nlat = CHARM(crd_point_get_local_nlat)(pnt);
    size_t local_nlon = CHARM(crd_point_get_local_nlon)(pnt);
    size_t npnt = local_nlat;


    if (local_nlat != local_nlon)
    {
        snprintf(err_msg, CHARM_ERR_MAX_MSG,
                         "The number of latitudes \"%zu\" do not "
                         "match the number of longitudes \"%zu\".",
                         local_nlat, local_nlon);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto BARRIER_1;
    }
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






    /* Synthesis */
    /* --------------------------------------------------------------------- */
    /* Get some constants */
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
        goto BARRIER_1;
    }


    if ((dr == GRAD_1) && (dlat == GRAD_1) && (dlon == GRAD_1))
        grad = 1;
    else if ((dr == GRAD_2) && (dlat == GRAD_2) && (dlon == GRAD_2))
        grad = 2;
    else
        grad = 0;


    /* ................................................................. */
    shcs_block = CHARM(shc_block_init)(shcs);
    CHECK_NULL(shcs_block, BARRIER_1);


    ips = (INT *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                       nmax * SIMD_SIZE * BLOCK_S,
                                       sizeof(INT));
    CHECK_NULL(ips, BARRIER_1);


    ps = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                       nmax * SIMD_SIZE * BLOCK_S,
                                       sizeof(REAL));
    CHECK_NULL(ps, BARRIER_1);


    tv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                       sizeof(REAL));
    CHECK_NULL(tv, BARRIER_1);


    uv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                       sizeof(REAL));
    CHECK_NULL(uv, BARRIER_1);


    lonv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                         SIMD_SIZE * BLOCK_S, sizeof(REAL));
    CHECK_NULL(lonv, BARRIER_1);


    pnt_rv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                           sizeof(REAL));
    CHECK_NULL(pnt_rv, BARRIER_1);


    tmpv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                         sizeof(REAL));
    CHECK_NULL(tmpv, BARRIER_1);


    size_t nrpows = BLOCK_S * (nmax + 2 + dorder);
    rpows = (REAL_SIMD *)CHARM(calloc_aligned)(SIMD_MEMALIGN, nrpows,
                                               sizeof(REAL_SIMD));
    CHECK_NULL(rpows, BARRIER_1);
    for (size_t i = 0; i < nrpows; i++)
        rpows[i] = SET1_R(PREC(1.0));
    /* ................................................................. */


    size_t ipv, l, idx;
    int err_glob = 0;
    unsigned long lc_err_glob = 0;


    MISC_SD_CALLOC_REAL_SIMD_ERR(t, BLOCK_S, SIMD_BLOCK_S, err, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_ERR(u, BLOCK_S, SIMD_BLOCK_S, err, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_ERR(fi, SHS_MAX_NPAR * BLOCK_S,
                                 SHS_MAX_NPAR * SIMD_BLOCK_S, err, BARRIER_1);
    REAL_SIMD pnt_r;
    REAL_SIMD tmp = SET_ZERO_R;


BARRIER_1:
    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto FAILURE_1;


    const size_t imax  = CHARM(shs_get_imax)(npnt, BLOCK_S, pnt);
    const size_t istep = SIMD_SIZE * BLOCK_S;


    for (size_t i = 0; i < imax; i += istep)
    {
        /* ------------------------------------------------------------- */
        for (l = 0; l < BLOCK_S; l++)
        {
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                ipv = i + l * SIMD_SIZE + v;
                if (ipv < npnt)
                {
                    tv[v] = SIN(pnt->lat[ipv]);
                    uv[v] = COS(pnt->lat[ipv]);
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


            t[l]  = LOAD_R(&tv[0]);
            u[l]  = LOAD_R(&uv[0]);
            pnt_r = LOAD_R(&pnt_rv[0]);
            if (!r_eq_rref)
                CHARM(shs_rpows)(pnt_r, rref, nmax + 1 + dorder, BLOCK_S,
                                 rpows + l);


            /* Prepare arrays for sectorial Legendre functions */
            CHARM(leg_func_prepare)(uv, ps + l * SIMD_SIZE * nmax,
                                    ips + l * SIMD_SIZE * nmax, dm, nmax);
        }
        /* ------------------------------------------------------------- */


        /* ------------------------------------------------------------- */
        for (l = 0; l < SHS_MAX_NPAR * BLOCK_S; l++)
            fi[l] = SET_ZERO_R;
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
shared(nmax, err, pt, t, u, rpows, ri, r, ips, ps, dorder, dr, dlat, dlon) \
shared(grad, r_eq_rref, shcs, shcs_block, lonv, fi, npar) \
shared(err_glob, lc_err_glob) \
private(l, idx) MPI_VARS
#endif
        {
        /* ------------------------------------------------------------- */
        REAL lontmp, m_real;
        REAL_SIMD clonim, slonim;
        /* ------------------------------------------------------------- */


        /* ------------------------------------------------------------- */
        int err_priv = 0;
        unsigned long lc_err_priv = 0;


        CHARM(lc) *lc    = NULL;
        REAL *anm        = NULL;
        REAL *bnm        = NULL;
        REAL *enm        = NULL;
        REAL *clonimv    = NULL;
        REAL *slonimv    = NULL;
        MISC_SD_CALLOC_REAL_SIMD_INIT(zeros);
        MISC_SD_CALLOC_REAL_SIMD_INIT(fi_thread);
        /* ------------------------------------------------------------- */


        /* Fill "shcs_block" with coefficients starting at degree "0".  If
         * compiling without the MPI support, "shc_block_get_coeffs" needs to
         * be called only once, so no further calls in relation to the
         * order-dependent loop are necessary.  "CHARM(shc_block_get_coeffs)"
         * is a collective call if compiled with MPI support, so should be
         * executed before the "anm", etc allocations.  Otherwise, a new
         * barrier would be needed. */
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


        clonimv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        CHECK_NULL_OMP(clonimv, err_priv, BARRIER_2);


        slonimv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                                sizeof(REAL));
        CHECK_NULL_OMP(slonimv, err_priv, BARRIER_2);


        MISC_SD_CALLOC_REAL_SIMD_ERR(zeros, BLOCK_S, SIMD_BLOCK_S, err,
                                     BARRIER_2);
        MISC_SD_CALLOC_REAL_SIMD_ERR(fi_thread, SHS_MAX_NPAR * BLOCK_S,
                                     SHS_MAX_NPAR * SIMD_BLOCK_S, err,
                                     BARRIER_2);
        for (l = 0; l < BLOCK_S; l++)
            zeros[l] = SET_ZERO_R;
        for (l = 0; l < SHS_MAX_NPAR * BLOCK_S; l++)
            fi_thread[l] = SET_ZERO_R;
        /* ------------------------------------------------------------- */


        /* Loop over harmonic orders */
        /* ------------------------------------------------------------- */
        unsigned long m = shcs_block->mfirst;  /* At this point, this value
                                                * will always be zero, because
                                                * of calling
                                                * "shc_block_get_coeffs"
                                                * above. */


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


          /* OpenMP requires the presence of the initializer in the for loop
           * that follows, even though "m" is defined */
#if HAVE_OPENMP
#pragma omp for schedule(dynamic)
#endif
            for (m = mmin; m <= mmax; m++)
            {
                /* Apply polar optimization if asked to do so */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, &u[0],
                                                         BLOCK_S, pt))
                    continue;


                /* "anm" and "bnm" coefficients for Legendre recurrence
                 * relations and their derivatives */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);
                if (dorder > 0)
                    CHARM(leg_func_enm)(nmax, m, r, ri, enm);


                /* Computation of the lumped coefficients */
                /* --------------------------------------------------------- */
                /* Summation over harmonic degrees.  Note that the symmetry of
                 * Legendre functions cannot be utilized with scattered points,
                 * hence "&zeros[0]".  The "lc->a2" and "lc->b2" are not used
                 * with scattered points. */
#undef KERNEL_IO_PARS
#define KERNEL_IO_PARS (nmax, m, shcs_block, r_eq_rref, anm, bnm, enm,        \
                        &t[0], &u[0], ps, ips, &rpows[0], &zeros[0],          \
                        &zeros[0], dorder, lc)
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


                /* The longitudinal part of the synthesis */
                /* ......................................................... */
                m_real = (REAL)m;


                for (l = 0; l < BLOCK_S; l++)
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
                        idx += BLOCK_S;
                        LC_CS(r);


                        idx += BLOCK_S;
                        LC_CS(p);
                    }


                    if (grad > 1)
                    {
                        idx += BLOCK_S;
                        LC_CS(rr);


                        idx += BLOCK_S;
                        LC_CS(rp);


                        idx += BLOCK_S;
                        LC_CS(pp);
                    }
                }
                /* ......................................................... */


            } /* End of the loop over harmonic orders */


            /* If using OpenMP, we do not know now what is the value of "m",
             * because it likely differs across the OpenMP threads.  So we have
             * to manually set it, so that the next iteration of the "do" loop
             * can be properly processed. */
            m = mmax + 1;
        }
        while (m <= nmax);


#if HAVE_OPENMP
#pragma omp critical
#endif
        {
        for (l = 0; l < npar * BLOCK_S; l++)
            fi[l] = ADD_R(fi[l], fi_thread[l]);


        lc_err_glob += lc_err_priv;
        }


FAILURE_2:
        CHARM(shs_lc_free)(lc);
        free(anm);
        free(bnm);
        free(enm);
        CHARM(free_aligned)(clonimv);
        CHARM(free_aligned)(slonimv);
        MISC_SD_FREE(zeros);
        MISC_SD_FREE(fi_thread);
        }
        /* ------------------------------------------------------------- */


        /* Final part of the synthesis */
        for (size_t p = 0; p < npar; p++)
            CHARM(shs_sctr_mulc)(i, npnt, pnt->type, mur, tmp, tmpv,
                                 &fi[p * BLOCK_S], f[p]);


    } /* End of the loop over the evaluation points */


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
    CHARM(shc_block_free)(shcs_block);
    CHARM(free_aligned)(ips);
    CHARM(free_aligned)(ps);
    CHARM(free_aligned)(tv);
    CHARM(free_aligned)(uv);
    CHARM(free_aligned)(lonv);
    CHARM(free_aligned)(pnt_rv);
    CHARM(free_aligned)(tmpv);
    CHARM(free_aligned)(rpows);
    MISC_SD_FREE(t);
    MISC_SD_FREE(u);
    MISC_SD_FREE(fi);
    /* --------------------------------------------------------------------- */






    return;
}
