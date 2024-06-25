/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "../prec.h"
#include "../shc/shc_reset_coeffs.h"
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../leg/leg_func_xnum.h"
#include "../leg/leg_func_use_xnum.h"
#include "../crd/crd_grd_check_symm.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../misc/misc_is_nearly_equal.h"
#include "../misc/misc_polar_optimization_threshold.h"
#include "../misc/misc_polar_optimization_apply.h"
#if CHARM_OPENMP
#   include <omp.h>
#endif
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
#define LOOP_ITER(n, a, b)                                                    \
    anms = SET1_R(anm[(n)]);                                                  \
    bnms = SET1_R(bnm[(n)]);                                                  \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        PNM_RECURRENCE(x[l], y[l], pnm2[l], t[l], anms, bnms);                \
        RECURRENCE_NEXT_ITER(y[l], x[l], pnm2[l]);                            \
    }                                                                         \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        shcs->c[m][idx] += SUM_R(MUL_R(pnm2[l], a[l]));                       \
    }                                                                         \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        shcs->s[m][idx] += SUM_R(MUL_R(pnm2[l], b[l]));                       \
    }                                                                         \
                                                                              \
                                                                              \
    idx++;
/* ------------------------------------------------------------------------- */






void CHARM(sha_point)(const CHARM(point) *pnt, const REAL *f,
                      unsigned long nmax, CHARM(shc) *shcs, CHARM(err) *err)
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
    const int pnt_type = pnt->type;
    unsigned long nmax_grd;


    /* Get the maximum degree for which the grid in "pnt" was created.  The
     * value is derived from the number of latitudes "pnt_nlat". */
    if (pnt_type == CHARM_CRD_POINT_GRID_GL)
    {
        /* In case of the Gauss--Legendre grid, it holds that "nmax_grd
         * = pnt_nlat - 1". */
        nmax_grd = pnt_nlat - 1;
    }
    else if ((pnt_type == CHARM_CRD_POINT_GRID_DH1) ||
             (pnt_type == CHARM_CRD_POINT_GRID_DH2))
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
    const REAL r0 = pnt->r[0];
    for (size_t i = 1; i < pnt->nlat; i++)
    /* In the line above, we intentionally used "pnt->nlat" instead of
     * "pnt_nlat", since the latter may have been modified, depending on the
     * grid type. */
    {
        if (!CHARM(misc_is_nearly_equal)(pnt->r[i], r0, CHARM(glob_threshold)))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "All spherical radii in \"pnt->r\" must be "
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
    if (pnt_nlat % 2)
    {
        even = 0; /* The number of latitudes is an odd number. The grid could
                   * possibly contain the zero latitude */
        nlatdo = (pnt_nlat + 1) / 2;
    }
    else
    {
        even = 1; /* The number of latitudes is an even number. The grid does
                   * not contain the zero latitude */
        nlatdo = pnt_nlat / 2;
    }
    /* --------------------------------------------------------------------- */






    /* Get the number of longitudes in "pnt" */
    /* --------------------------------------------------------------------- */
    const size_t pnt_nlon     = pnt->nlon;
    const size_t pnt_nlon_fft = pnt_nlon / 2 + 1;
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    int FAILURE_glob = 0;
    REAL *r                  = NULL;
    REAL *ri                 = NULL;
    REAL *dm                 = NULL;
    REAL *ftmp_in            = NULL;
    FFTWC(complex) *ftmp_out = NULL;
    FFTW(plan) plan          = NULL;
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






    /* Auxiliary variable entering the computation of the lumped coefficients
     * */
    /* --------------------------------------------------------------------- */
    REAL c = PREC(0.0);


    /* Here, we have to use the maximum degree "nmax_grd" for which the input
     * grid in "pnt" was created. */
    if ((pnt_type == CHARM_CRD_POINT_GRID_GL) ||
        (pnt_type == CHARM_CRD_POINT_GRID_DH1))
        c = PI / (REAL)(nmax_grd + 1);
    else if (pnt_type == CHARM_CRD_POINT_GRID_DH2)
        c = PI / (REAL)(2 * nmax_grd + 2);
    /* --------------------------------------------------------------------- */






    /* Create a plan for FFT */
    /* --------------------------------------------------------------------- */
#if CHARM_OPENMP && FFTW3_OMP
    if (FFTW(init_threads)() == 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFFTWINIT,
                       "FFTW failed to initialize threads.");
        return;
    }


    FFTW(plan_with_nthreads)(omp_get_max_threads());
#endif
    ftmp_in = (REAL *)FFTW(malloc)(pnt_nlon * sizeof(REAL));
    if (ftmp_in == NULL)
    {
        FFTW(free)(ftmp_in);
        FAILURE_glob = 1;
        goto FAILURE;
    }
    ftmp_out = (FFTWC(complex) *)FFTW(malloc)(pnt_nlon_fft *
                                              sizeof(FFTWC(complex)));
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






    /* --------------------------------------------------------------------- */
    /* Set all coefficients in "shcs" to zero */
    CHARM(shc_reset_coeffs)(shcs);


    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


    {
        REAL_SIMD pnm0[SIMD_BLOCK], pnm1[SIMD_BLOCK], pnm2[SIMD_BLOCK];
        REAL_SIMD x[SIMD_BLOCK], y[SIMD_BLOCK], z[SIMD_BLOCK];
        REAL_SIMD t[SIMD_BLOCK], u[SIMD_BLOCK];
        RI_SIMD   ix[SIMD_BLOCK], iy[SIMD_BLOCK], iz[SIMD_BLOCK];
        RI_SIMD   ixy[SIMD_BLOCK];
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
        ABS_R_INIT;
#endif
        REAL_SIMD symm[SIMD_BLOCK], latsin[SIMD_BLOCK];
        REAL_SIMD am[SIMD_BLOCK],   bm[SIMD_BLOCK];
        REAL_SIMD a2m[SIMD_BLOCK],  b2m[SIMD_BLOCK];


        /* ................................................................. */
        int *ips  = NULL;
        REAL *ps  = NULL;
#if !(CHARM_OPENMP)
        REAL *anm = NULL;
        REAL *bnm = NULL;
#endif
        REAL *tv      = NULL;
        REAL *uv      = NULL;
        REAL *symmv   = NULL;  /* To be store only "1.0" or "0.0" */
        REAL *latsinv = NULL;  /* To be store only "1.0" or "0.0" */
        REAL *a       = NULL;
        REAL *b       = NULL;
        REAL *a2      = NULL;
        REAL *b2      = NULL;
        REAL *ftmp_in = NULL;
        FFTWC(complex) *ftmp_out = NULL;


        ips = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           SIMD_SIZE * SIMD_BLOCK * nmax,
                                           sizeof(int));
        if (ips == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ps = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           SIMD_SIZE * SIMD_BLOCK * nmax,
                                           sizeof(REAL));
        if (ps == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
#if !(CHARM_OPENMP)
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
        tv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                           sizeof(REAL));
        if (tv == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        uv = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                           sizeof(REAL));
        if (uv == NULL)
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
                                          SIMD_SIZE * SIMD_BLOCK *
                                          pnt_nlon_fft, sizeof(REAL));
        if (a == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                          SIMD_SIZE * SIMD_BLOCK *
                                          pnt_nlon_fft, sizeof(REAL));
        if (b == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        a2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           SIMD_SIZE * SIMD_BLOCK *
                                           pnt_nlon_fft, sizeof(REAL));
        if (a2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           SIMD_SIZE * SIMD_BLOCK *
                                           pnt_nlon_fft, sizeof(REAL));
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
        ftmp_out = (FFTWC(complex) *)FFTW(malloc)(pnt_nlon_fft *
                                                  sizeof(FFTWC(complex)));
        if (ftmp_out == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        /* ................................................................. */



        REAL_SIMD anms, bnms;
        REAL_SIMD amp[SIMD_BLOCK], amm[SIMD_BLOCK];
        REAL_SIMD bmp[SIMD_BLOCK], bmm[SIMD_BLOCK];
        REAL cw;
        REAL_SIMD wlf;
        _Bool npm_even; /* True if "n + m" is even */
        size_t l, ipv; /* "i + v" */
        _Bool ds[SIMD_BLOCK]; /* Dynamical switching */


        /* Loop over latitudes */
        for (size_t i = 0; i < SIMD_MULTIPLE(nlatdo, SIMD_SIZE * SIMD_BLOCK);
             i += SIMD_SIZE * SIMD_BLOCK)
        {
            for (l = 0; l < SIMD_BLOCK; l++)
            {
                for (size_t v = 0; v < SIMD_SIZE; v++)
                {
                    /* Check whether the symmetry property of LFs needs to be
                     * applied */
                    /* ----------------------------------------------------- */
                    ipv = i + l * SIMD_SIZE + v;
                    CHARM(crd_grd_check_symm)(ipv, v, pnt_type, nlatdo, 0,
                                              even, symmv, latsinv);
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
                    memcpy(ftmp_in, f + ipv * pnt_nlon,
                           pnt_nlon * sizeof(REAL));
                    FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


                    cw = c * pnt->w[ipv];
                    for (size_t j = 0; j < pnt_nlon_fft; j++)
                    {
                        a[j * (SIMD_SIZE * SIMD_BLOCK) +
                          l * SIMD_SIZE + v]  =  cw * ftmp_out[j][0];
                        b[j * (SIMD_SIZE * SIMD_BLOCK) +
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
                            a2[j * (SIMD_SIZE * SIMD_BLOCK) +
                               l * SIMD_SIZE + v]  =  cw * ftmp_out[j][0];
                            b2[j * (SIMD_SIZE * SIMD_BLOCK) +
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
#if CHARM_OPENMP


#   undef SIMD_VARS
#   ifdef SIMD
#       define SIMD_VARS shared(zero_ri, one_ri, mone_ri, zero_r) \
                         shared(BIG_r, BIGI_r, BIGS_r, BIGSI_r) \
                         shared(NONSIGNBITS_R) \
                         private(tmp1_r, tmp2_r, mask1, mask2, mask3)
#   else
#       define SIMD_VARS
#   endif


#pragma omp parallel default(none) \
shared(nmax, symm, r, ri, a, b, a2, b2, shcs, t, u, ps, ips) \
shared(latsin, pt, ROOT3_r, FAILURE_glob, err) \
private(am, bm, a2m, b2m, amp, amm, bmp, bmm, anms, bnms) \
private(x, ix, y, iy, wlf, ixy, z, iz) \
private(pnm0, pnm1, pnm2, npm_even, ds, l) SIMD_VARS
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
                if (CHARM(err_isempty)(err))
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EMEM, CHARM_ERR_MALLOC_FAILURE);


                /* OK, and now all threads go to the "FAILURE_2_parallel"
                 * label to deallocate all the memory that might be allocated
                 * before the allocation failure. */
                goto FAILURE_2_parallel;
            }
            /* ............................................................. */


            /* Loop over harmonic orders */
            unsigned long m;
#pragma omp for schedule(dynamic) private(m)
#endif
            for (m = 0; m <= nmax; m++)
            {

                /* Apply polar optimization if asked to do so */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, &u[0],
                                                         SIMD_BLOCK, pt))
                    continue;


                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Some useful substitutions */
                /* --------------------------------------------------------- */
                for (l = 0; l < SIMD_BLOCK; l++)
                {
                    am[l]  = LOAD_R(&a[(SIMD_SIZE * SIMD_BLOCK) * m +
                                       l * SIMD_SIZE]);
                    bm[l]  = LOAD_R(&b[(SIMD_SIZE * SIMD_BLOCK) * m +
                                       l * SIMD_SIZE]);
                    a2m[l] = LOAD_R(&a2[(SIMD_SIZE * SIMD_BLOCK) * m +
                                        l * SIMD_SIZE]);
                    b2m[l] = LOAD_R(&b2[(SIMD_SIZE * SIMD_BLOCK) * m +
                                        l * SIMD_SIZE]);


                    amp[l] = MUL_R(ADD_R(am[l], MUL_R(symm[l], a2m[l])),
                                   latsin[l]);
                    amm[l] = MUL_R(SUB_R(am[l], MUL_R(symm[l], a2m[l])),
                                   latsin[l]);
                    bmp[l] = MUL_R(ADD_R(bm[l], MUL_R(symm[l], b2m[l])),
                                   latsin[l]);
                    bmm[l] = MUL_R(SUB_R(bm[l], MUL_R(symm[l], b2m[l])),
                                   latsin[l]);
                }
                /* --------------------------------------------------------- */


                /* Computation of spherical harmonic coefficients */
                if (m == 0)
                {
                    /* Zonal harmonics */
                    /* ----------------------------------------------------- */

                    /* P00 */
                    for (l = 0; l < SIMD_BLOCK; l++)
                    {
                        pnm0[l] = SET1_R(PREC(1.0));
                        /* C00 */
                        shcs->c[0][0] += SUM_R(MUL_R(pnm0[l], amp[l]));
                    }


                    /* P10 */
                    if (nmax >= 1)
                    {
                        for (l = 0; l < SIMD_BLOCK; l++)
                        {
                            pnm1[l] = MUL_R(ROOT3_r, t[l]);
                            /* C10 */
                            shcs->c[0][1] += SUM_R(MUL_R(pnm1[l], amm[l]));
                        }
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


                            for (l = 0; l < SIMD_BLOCK; l++)
                            {
                                pnm2[l] = SUB_R(MUL_R(MUL_R(anms,
                                                            t[l]), pnm1[l]),
                                                MUL_R(bnms,
                                                      pnm0[l]));
                                pnm0[l] = pnm1[l];
                                pnm1[l] = pnm2[l];
                            }


                            for (l = 0; l < SIMD_BLOCK; l++)
                            {
                                /* C20, C30, ..., Cnmax,0 */
                                shcs->c[0][n] += SUM_R(MUL_R(pnm2[l],
                                                       (npm_even) ? amp[l] :
                                                                    amm[l]));
                            }
                        }
                    }
                    /* ----------------------------------------------------- */

                }
                else /* Non-zonal harmonics */
                {

                    /* Sectorial harmonics */
                    /* ----------------------------------------------------- */
                    for (l = 0; l < SIMD_BLOCK; l++)
                    {
#ifdef SIMD
                        PNM_SECTORIAL_XNUM_SIMD(x[l], ix[l],
                                                ps[(SIMD_SIZE * nmax) * l +
                                                   (m - 1) * SIMD_SIZE],
                                                ips[(SIMD_SIZE * nmax) * l +
                                                   (m - 1) * SIMD_SIZE],
                                                pnm0[l], BIG_r,
                                                zero_r, zero_ri,
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
                    for (l = 0; l < SIMD_BLOCK; l++)
                    {
                        shcs->c[m][0] += SUM_R(MUL_R(pnm0[l], amp[l]));
                        shcs->s[m][0] += SUM_R(MUL_R(pnm0[l], bmp[l]));
                    }
                    /* ----------------------------------------------------- */


                    /* Tesseral harmonics */
                    /* ----------------------------------------------------- */
                    if (m < nmax)
                    {
                        anms = SET1_R(anm[m + 1]);
                        bnms = SET1_R(bnm[m + 1]);


                        for (l = 0; l < SIMD_BLOCK; l++)
                        {
#ifdef SIMD
                            PNM_SEMISECTORIAL_XNUM_SIMD(x[l], y[l],
                                                        ix[l], iy[l],
                                                        wlf, t[l],
                                                        anms, pnm1[l],
                                                        mask1, mask2, mask3,
                                                        zero_r, zero_ri,
                                                        mone_ri,
                                                        BIG_r, BIGS_r,  BIGI_r,
                                                        SEMISECTORIALS);
#else
                            PNM_SEMISECTORIAL_XNUM(x[l], y[l], ix[l], iy[l],
                                                   wlf, t[l], anms,
                                                   pnm1[l]);
#endif
                        }


                        /* Cm+1,m; Sm+1,m */
                        for (l = 0; l < SIMD_BLOCK; l++)
                        {
                            shcs->c[m][1] += SUM_R(MUL_R(pnm1[l], amm[l]));
                            shcs->s[m][1] += SUM_R(MUL_R(pnm1[l], bmm[l]));
                        }


                        /* Loop over degrees */
                        /* ------------------------------------------------- */
                        for (l = 0; l < SIMD_BLOCK; l++)
                            ds[l] = 0;


                        /* Is "n + m" even?  Since we start the loop with "n
                         * = m + 2", then the parity of the first "m + 2 + m"
                         * is always even.  Then, it changes with every loop
                         * iteration. */
                        npm_even = 1;


                        unsigned long n;
                        unsigned long idx = 2;
                        for (n = (m + 2);
                             CHARM(leg_func_use_xnum(ds, SIMD_BLOCK)) &&
                             n <= nmax;
                             n++, npm_even = !npm_even, idx++)
                        {
                            anms = SET1_R(anm[n]);
                            bnms = SET1_R(bnm[n]);


                            /* Compute tesseral Legendre function */
                            for (l = 0; l < SIMD_BLOCK; l++)
                            {
#ifdef SIMD
                                PNM_TESSERAL_XNUM_SIMD(x[l], y[l], z[l],
                                                       ix[l], iy[l], iz[l],
                                                       ixy[l],
                                                       wlf, t[l],
                                                       anms, bnms,
                                                       pnm2[l],
                                                       tmp1_r, tmp2_r,
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
                            if (npm_even)
                            {
                                for (l = 0; l < SIMD_BLOCK; l++)
                                     shcs->c[m][idx] += SUM_R(MUL_R(pnm2[l],
                                                                    amp[l]));
                                for (l = 0; l < SIMD_BLOCK; l++)
                                     shcs->s[m][idx] += SUM_R(MUL_R(pnm2[l],
                                                                    bmp[l]));
                            }
                            else
                            {
                                for (l = 0; l < SIMD_BLOCK; l++)
                                     shcs->c[m][idx] += SUM_R(MUL_R(pnm2[l],
                                                                    amm[l]));
                                for (l = 0; l < SIMD_BLOCK; l++)
                                     shcs->s[m][idx] += SUM_R(MUL_R(pnm2[l],
                                                                    bmm[l]));
                            }


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


#if CHARM_OPENMP
FAILURE_2_parallel:
            free(anm); free(bnm);
            }
#endif
            /* ------------------------------------------------------------- */


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


FAILURE_1:
        if ((FAILURE_glob != 0) && CHARM(err_isempty(err)))
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);


        FFTW(free)(ftmp_in);        FFTW(free)(ftmp_out);
        CHARM(free_aligned)(tv);    CHARM(free_aligned)(uv);
        CHARM(free_aligned)(symmv); CHARM(free_aligned)(latsinv);
        CHARM(free_aligned)(ips);   CHARM(free_aligned)(ps);
        CHARM(free_aligned)(a);     CHARM(free_aligned)(b);
        CHARM(free_aligned)(a2);    CHARM(free_aligned)(b2);
#if !(CHARM_OPENMP)
        free(anm); free(bnm);
#endif


    }
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if ((FAILURE_glob != 0) && CHARM(err_isempty(err)))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    free(r); free(ri); free(dm);
    FFTW(destroy_plan)(plan);
#if CHARM_OPENMP && FFTW3_OMP
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
    const REAL c2 = PREC(1.0) / (PREC(4.0) * PI) * (r0 / shcs->mu);
    unsigned long m;
#if CHARM_OPENMP
    #pragma omp parallel for default(none) shared(shcs, nmax, c2) private(m)
#endif
    for (m = 0; m <= nmax; m++)
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
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    return;
}

