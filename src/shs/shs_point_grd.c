/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "../prec.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_enm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../crd/crd_grd_check_symm.h"
#include "../crd/crd_point_isQuadGrid.h"
#include "../crd/crd_point_isCustGrid.h"
#include "../misc/misc_arr_chck_lin_incr.h"
#include "../misc/misc_arr_chck_symm.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
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
#include "shs_max_npar.h"
#include "shs_point_grd.h"
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


    /* Finally, if the grid is symmetric, we modify the number of latitudes
     * "pnt_nlat" to be equal to the number of latitudes on one hemisphere
     * only (including the equator if present). This is because the
     * time-consuming "for loop" over evaluation points now needs to run for
     * one hemisphere only, while the results for the other hemisphere are
     * obtained by exploiting the symmetry property of Legendre functions. This
     * reduces the number of Legendre functions that need to be evaluated by
     * a factor of ~2, so saves some computational time */
    size_t nlatdo = (symm) ? (pnt_nlat + 1 - even) / 2 : pnt_nlat;
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
    REAL deltalon = (pnt_nlon > 1) ? pnt->lon[1] - pnt->lon[0] : PREC(0.0);


    /* Auxiliary constant to be used only in case the PSLR algorithm is
     * applied.  To suppress, a compiler warning, the value is initialized to
     * zero. */
    REAL lon0 = PREC(0.0);


    /* Length of the lumped coefficients arrays in case FFT will be applied */
    size_t nfc = pnt_nlon / 2 + 1;


    _Bool use_fft = CHARM(shs_grd_point_fft_check)(pnt, deltalon, nmax);
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
    _Bool r_eq_rref = CHARM(shs_r_eq_rref)(pnt, shcs);
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    int FAILURE_glob       = 0;
    REAL *r                = NULL;
    REAL *ri               = NULL;
    REAL *dm               = NULL;
    FFTWC(complex) *fc_tmp = NULL;
    REAL *ftmp             = NULL;
    FFTW(plan) plan        = NULL;
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
    /* --------------------------------------------------------------------- */






    /* Create a FFT plan */
    /* --------------------------------------------------------------------- */
    if (use_fft)
    {
        fc_tmp = (FFTWC(complex) *)FFTW(malloc)(nfc * sizeof(FFTWC(complex)));
        if (fc_tmp == NULL)
        {
            FFTW(free)(fc_tmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }
        ftmp = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
        if (ftmp == NULL)
        {
            FFTW(free)(fc_tmp);
            FFTW(free)(ftmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }


        plan = FFTW(plan_dft_c2r_1d)(pnt_nlon, fc_tmp, ftmp, FFTW_ESTIMATE);
        if (plan == NULL)
        {
            FFTW(free)(fc_tmp);
            FFTW(free)(ftmp);
            FAILURE_glob = 1;
            goto FAILURE;
        }


        FFTW(free)(fc_tmp); FFTW(free)(ftmp);
    }
    /* --------------------------------------------------------------------- */






    /* Loop over grid latitudes */
    /* --------------------------------------------------------------------- */
    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


    /* Radius of the reference sphere that is associated with the spherical
     * harmonic coefficients */
    REAL_SIMD rref = SET1_R(shcs->r);


#if CHARM_OPENMP
#pragma omp parallel default(none) \
shared(f, shcs, nmax, pnt, pnt_nlat, pnt_nlon, dm, r, ri, nlatdo, pnt_type) \
shared(lon0, deltalon, even, symm, FAILURE_glob, mur, err, nfc, plan) \
shared(use_fft, pt, rref, r_eq_rref, dr, dlat, dlon, dorder, npar, grad)
#endif
    {
        /* ................................................................. */
        /* An indicator for failed memory initializations on each thread,
         * a private variable. */
        int FAILURE_priv = 0;


        INT *ips            = NULL;
        REAL *ps            = NULL;
        REAL *latv          = NULL;
        REAL *tv            = NULL;
        REAL *uv            = NULL;
        REAL *symmv         = NULL;
        REAL *latsinv       = NULL;
        REAL *pnt_rv        = NULL;
        REAL *pnt_r2v       = NULL;
        REAL *fc_simd       = NULL;
        REAL *fc2_simd      = NULL;
        REAL *anm           = NULL;
        REAL *bnm           = NULL;
        REAL *enm           = NULL;
        FFTWC(complex) *fc  = NULL;
        FFTWC(complex) *fc2 = NULL;
        REAL *ftmp          = NULL;
        REAL *fi            = NULL;
        REAL *fi2           = NULL;


        size_t nfi_1par = pnt_nlon * SIMD_SIZE * SIMD_BLOCK_S;
        size_t nfi = npar * nfi_1par;


        ips = (INT *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           nmax * SIMD_SIZE * SIMD_BLOCK_S,
                                           sizeof(INT));
        if (ips == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        ps = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           nmax * SIMD_SIZE * SIMD_BLOCK_S,
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
        symmv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                              SIMD_SIZE * SIMD_BLOCK_S,
                                              sizeof(REAL));
        if (symmv == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }
        latsinv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                SIMD_SIZE * SIMD_BLOCK_S,
                                                sizeof(REAL));
        if (latsinv == NULL)
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


        if (use_fft)
        {
            ftmp = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
            if (ftmp == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
            fc = (FFTWC(complex) *)FFTW(malloc)(nfc * sizeof(FFTWC(complex)));
            if (fc == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
            memset(fc, 0, nfc * sizeof(FFTWC(complex)));
            fc_simd = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                    npar * nfc * SIMD_SIZE *
                                                    SIMD_BLOCK_S * 2,
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


        if (symm)
        {
            pnt_r2v = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                    SIMD_SIZE, sizeof(REAL));
            if (pnt_r2v == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }


            if (use_fft)
            {
                fc2 = (FFTWC(complex) *)FFTW(malloc)(nfc *
                                                     sizeof(FFTWC(complex)));
                if (fc2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
                memset(fc2, 0, nfc * sizeof(FFTWC(complex)));
                fc2_simd = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                         npar * nfc *
                                                         SIMD_SIZE *
                                                         SIMD_BLOCK_S * 2,
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


        size_t ipv, l;
        const size_t size_blk2 = SIMD_SIZE * SIMD_BLOCK_S * 2;


        REAL_SIMD t[SIMD_BLOCK_S], u[SIMD_BLOCK_S], symm_simd[SIMD_BLOCK_S];
        REAL_SIMD pnt_r[SIMD_BLOCK_S], pnt_r2[SIMD_BLOCK_S];
        REAL_SIMD ratio[SIMD_BLOCK_S], ratio2[SIMD_BLOCK_S];
        REAL_SIMD ratiom[SIMD_BLOCK_S], ratio2m[SIMD_BLOCK_S];
        for (l = 0; l < SIMD_BLOCK_S; l++)
            ratio[l] = ratio2[l] = ratiom[l] = ratio2m[l] = SET_ZERO_R;
        CHARM(lc) lc;
        CHARM(shs_lc_init)(&lc);


        size_t i;
#if CHARM_OPENMP
#pragma omp for schedule(dynamic) private(i)
#endif
        for (i = 0; i < SIMD_MULTIPLE(nlatdo, SIMD_SIZE * SIMD_BLOCK_S);
             i += SIMD_SIZE * SIMD_BLOCK_S)
        {
            for (l = 0; l < SIMD_BLOCK_S; l++)
            {
                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    /* Check whether the symmetry property of LFs can be
                     * applied */
                    /* ----------------------------------------------------- */
                    ipv = i + l * SIMD_SIZE + v;
                    CHARM(crd_grd_check_symm)(ipv, v, pnt_type, nlatdo, symm,
                                              even, symmv + l * SIMD_SIZE,
                                              latsinv + l * SIMD_SIZE);


                    if (latsinv[l * SIMD_SIZE + v] == 1)
                    {
                        latv[v]   = pnt->lat[ipv];
                        tv[v]     = SIN(pnt->lat[ipv]);
                        uv[v]     = COS(pnt->lat[ipv]);
                        pnt_rv[v] = pnt->r[ipv];


                        if (symm)
                        {
                            if (((pnt_type == CHARM_CRD_POINT_GRID_DH1) ||
                                 (pnt_type == CHARM_CRD_POINT_GRID_DH2)) &&
                                (ipv == 0))
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
                pnt_r[l]     = LOAD_R(&pnt_rv[0]);
                symm_simd[l] = LOAD_R(&symmv[l * SIMD_SIZE]);


                ratio[l]  = DIV_R(rref, pnt_r[l]);
                ratiom[l] = ratio[l];
                if (symm)
                {
                    pnt_r2[l]  = LOAD_R(&pnt_r2v[0]);
                    ratio2[l]  = DIV_R(rref, pnt_r2[l]);
                    ratio2m[l] = ratio2[l];
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


            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Apply polar optimization if asked to do so */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, &u[0],
                                                         SIMD_BLOCK_S, pt))
                    goto UPDATE_RATIOS;


                /* "anm" and "bnm" coefficients for Legendre recurrence
                 * relations and for their derivatives */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);
                if (dorder > 0)
                    CHARM(leg_func_enm)(nmax, m, r, ri, enm);



                /* Computation of the lumped coefficients */
                /* --------------------------------------------------------- */
#undef KERNEL_IO_PARS
#define KERNEL_IO_PARS (nmax, m, shcs, r_eq_rref, anm, bnm, enm,              \
                        &t[0], &u[0], ps, ips,                                \
                        &ratio[0], &ratio2[0], &ratiom[0], &ratio2m[0],       \
                        &symm_simd[0], dorder, &lc);
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


                /* --------------------------------------------------------- */
                if (use_fft)
                    CHARM(shs_grd_fft_lc)(m, deltalon, grad, &lc,
                                          symm, &symm_simd[0], pnt_type,
                                          nfc, fc_simd, fc2_simd);
                else
                    CHARM(shs_grd_lr)(m, lon0, deltalon, pnt_nlon, pnt_type,
                                      grad, nfi_1par, &lc, symm,
                                      fi, fi2);
                /* --------------------------------------------------------- */


UPDATE_RATIOS:
                for (l = 0; l < SIMD_BLOCK_S; l++)
                {
                    ratiom[l] = MUL_R(ratiom[l], ratio[l]);
                    if (symm)
                        ratio2m[l] = MUL_R(ratio2m[l], ratio2[l]);
                }


            } /* End of the loop over harmonic orders */
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
        /* ----------------------------------------------------------------- */


        /* Free the heap memory */
        /* ----------------------------------------------------------------- */
FAILURE_2_parallel:
        CHARM(free_aligned)(ips);       CHARM(free_aligned)(ps);
        CHARM(free_aligned)(latv);
        CHARM(free_aligned)(tv);        CHARM(free_aligned)(uv);
        CHARM(free_aligned)(symmv);     CHARM(free_aligned)(latsinv);
        CHARM(free_aligned)(pnt_rv);    CHARM(free_aligned)(pnt_r2v);
        CHARM(free_aligned)(fi);        CHARM(free_aligned)(fi2);
        CHARM(free_aligned)(fc_simd);   CHARM(free_aligned)(fc2_simd);
        free(anm);
        free(bnm);
        free(enm);
        FFTW(free)(fc); FFTW(free)(fc2);
        FFTW(free)(ftmp);
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
