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
    int pnt_type = pnt->type;
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
    REAL r0 = pnt->r[0];
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
    size_t pnt_nlon     = pnt->nlon;
    size_t pnt_nlon_fft = pnt_nlon / 2 + 1;
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
    ftmp_out = (FFTW(complex) *)FFTW(malloc)(pnt_nlon_fft *
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






    /* --------------------------------------------------------------------- */
    /* Set all coefficients in "shcs" to zero */
    CHARM(shc_reset_coeffs)(shcs);


    /* Get the polar optimization threshold */
    REAL_SIMD pt = CHARM(misc_polar_optimization_threshold)(nmax);


    {
        REAL_SIMD pnm0, pnm1, pnm2;
        REAL_SIMD x, y, z, t, u;
        RI_SIMD   ix, iy, iz, ixy;
        REAL_SIMD ROOT3_r = SET1_R(ROOT3);
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
        REAL_SIMD symm, latsin;
        REAL_SIMD am, bm, a2m, b2m;


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
        FFTW(complex) *ftmp_out = NULL;


        ips = (int *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE * nmax,
                                           sizeof(int));
        if (ips == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ps = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE * nmax,
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
                                          SIMD_SIZE * pnt_nlon_fft,
                                          sizeof(REAL));
        if (a == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                          SIMD_SIZE * pnt_nlon_fft,
                                          sizeof(REAL));
        if (b == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        a2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           SIMD_SIZE * pnt_nlon_fft,
                                           sizeof(REAL));
        if (a2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b2 = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                           SIMD_SIZE * pnt_nlon_fft,
                                           sizeof(REAL));
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
        ftmp_out = (FFTW(complex) *)FFTW(malloc)(pnt_nlon_fft *
                                                 sizeof(FFTW(complex)));
        if (ftmp_out == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        /* ................................................................. */



        REAL_SIMD amp, amm, bmp, bmm;
        REAL cw;
        REAL_SIMD wlf;
        _Bool npm_even; /* True if "n + m" is even */
        size_t ipv; /* "i + v" */
        unsigned long nmm; /* "n - m" */
        _Bool ds; /* Dynamical switching */


        /* Loop over latitudes */
        for (size_t i = 0; i < SIMD_GET_MULTIPLE(nlatdo); i += SIMD_SIZE)
        {
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                /* Check whether the symmetry property of LFs needs to be
                 * applied */
                /* --------------------------------------------------------- */
                ipv = i + v;
                CHARM(crd_grd_check_symm)(ipv, v, pnt_type, nlatdo, 0, even,
                                          symmv, latsinv);
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
                /* --------------------------------------------------------- */


                /* Lumped coefficients for the southern hemisphere (including
                 * the equator) */
                /* --------------------------------------------------------- */
                memcpy(ftmp_in, f + ipv * pnt_nlon, pnt_nlon * sizeof(REAL));
                FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


                cw = c * pnt->w[ipv];
                for (size_t j = 0; j < pnt_nlon_fft; j++)
                {
                    a[j * SIMD_SIZE + v]  =  cw * ftmp_out[j][0];
                    b[j * SIMD_SIZE + v]  = -cw * ftmp_out[j][1];
                }
                /* --------------------------------------------------------- */


                /* Lumped coefficients for the northern hemisphere */
                /* --------------------------------------------------------- */
                if (symmv[v])
                {
                    memcpy(ftmp_in, f + (pnt_nlat - ipv - 1) * pnt_nlon,
                           pnt_nlon * sizeof(REAL));
                    FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


                    cw = c * pnt->w[pnt_nlat - ipv - 1];
                    for (size_t j = 0; j < pnt_nlon_fft; j++)
                    {
                        a2[j * SIMD_SIZE + v]  =  cw * ftmp_out[j][0];
                        b2[j * SIMD_SIZE + v]  = -cw * ftmp_out[j][1];
                    }
                }
                /* --------------------------------------------------------- */
            }


            t      = LOAD_R(&tv[0]);
            u      = LOAD_R(&uv[0]);
            symm   = LOAD_R(&symmv[0]);
            latsin = LOAD_R(&latsinv[0]);


            /* Prepare arrays for sectorial Legendre functions */
            /* ------------------------------------------------------------- */
            CHARM(leg_func_prepare)(uv, ps, ips, dm, nmax);
            /* ------------------------------------------------------------- */


            /* ------------------------------------------------------------- */
#if CHARM_OPENMP
#   ifdef SIMD
        #pragma omp parallel default(none) \
            shared(nmax, symm, r, ri, a, b, a2, b2, shcs, t, u, ps, ips) \
            shared(latsin, pt, ROOT3_r, FAILURE_glob, err) \
            private(am, bm, a2m, b2m, amp, amm, bmp, bmm) \
            private(x, ix, y, iy, wlf, ixy, z, iz) \
            private(pnm0, pnm1, pnm2, npm_even, nmm, ds) \
            shared(zero_ri, one_ri, mone_ri, zero_r) \
            shared(BIG_r, BIGI_r, BIGS_r, BIGSI_r, NONSIGNBITS_R) \
            private(tmp1_r, tmp2_r, mask1, mask2, mask3)
#   else
        #pragma omp parallel default(none) \
            shared(nmax, symm, r, ri, a, b, a2, b2, shcs, t, u, ps, ips) \
            shared(latsin, pt, ROOT3_r, FAILURE_glob, err) \
            private(am, bm, a2m, b2m, amp, amm, bmp, bmm) \
            private(x, ix, y, iy, wlf, ixy, z, iz) \
            private(pnm0, pnm1, pnm2, npm_even, nmm, ds)
#   endif
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
#pragma omp for schedule(dynamic)
#endif
            for (unsigned long m = 0; m <= nmax; m++)
            {

                /* Apply polar optimization if asked to do so */
                if (CHARM(misc_polar_optimization_apply)(m, nmax, u, pt))
                    continue;


                /* Computation of "anm" and "bnm" coefficients for Legendre
                 * recurrence relations */
                CHARM(leg_func_anm_bnm)(nmax, m, r, ri, anm, bnm);


                /* Some useful substitutions */
                /* --------------------------------------------------------- */
                am  = LOAD_R(&a[SIMD_SIZE * m]);
                bm  = LOAD_R(&b[SIMD_SIZE * m]);
                a2m = LOAD_R(&a2[SIMD_SIZE * m]);
                b2m = LOAD_R(&b2[SIMD_SIZE * m]);


                amp = MUL_R(ADD_R(am, MUL_R(symm, a2m)), latsin);
                amm = MUL_R(SUB_R(am, MUL_R(symm, a2m)), latsin);
                bmp = MUL_R(ADD_R(bm, MUL_R(symm, b2m)), latsin);
                bmm = MUL_R(SUB_R(bm, MUL_R(symm, b2m)), latsin);
                /* --------------------------------------------------------- */


                /* Computation of spherical harmonic coefficients */
                if (m == 0)
                {
                    /* Zonal harmonics */
                    /* ----------------------------------------------------- */

                    /* P00 */
                    pnm0 = SET1_R(PREC(1.0));
                    /* C00 */
                    shcs->c[0][0] += SUM_R(MUL_R(pnm0, amp));


                    /* P10 */
                    if (nmax >= 1)
                    {
                        pnm1 = MUL_R(ROOT3_r, t);
                        /* C10 */
                        shcs->c[0][1] += SUM_R(MUL_R(pnm1, amm));
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
                            pnm2 = SUB_R(MUL_R(MUL_R(SET1_R(anm[n]), t), pnm1),
                                         MUL_R(SET1_R(bnm[n]), pnm0));
                            /* C20, C30, ..., Cnmax,0 */
                            if (npm_even)
                                shcs->c[0][n] += SUM_R(MUL_R(pnm2, amp));
                            else
                                shcs->c[0][n] += SUM_R(MUL_R(pnm2, amm));


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
#ifdef SIMD
                    PNM_SECTORIAL_XNUM_SIMD(x, ix,
                                            ps[SIMD_SIZE * (m - 1)],
                                            ips[SIMD_SIZE * (m - 1)],
                                            pnm0, BIG_r, zero_r, zero_ri,
                                            mone_ri, mask1, mask2, SECTORIALS);
#else
                    PNM_SECTORIAL_XNUM(x, ix, ps[m - 1], ips[m - 1], pnm0);
#endif


                    /* Cm,m; Sm,m */
                    shcs->c[m][0] += SUM_R(MUL_R(pnm0, amp));
                    shcs->s[m][0] += SUM_R(MUL_R(pnm0, bmp));
                    /* ----------------------------------------------------- */


                    /* Tesseral harmonics */
                    /* ----------------------------------------------------- */
                    if (m < nmax)
                    {
#ifdef SIMD
                        PNM_SEMISECTORIAL_XNUM_SIMD(x, y, ix, iy, wlf, t,
                                                    anm[m + 1], pnm1,
                                                    mask1, mask2, mask3,
                                                    zero_r, zero_ri, mone_ri,
                                                    BIG_r, BIGS_r,  BIGI_r,
                                                    SEMISECTORIALS);
#else
                        PNM_SEMISECTORIAL_XNUM(x, y, ix, iy, wlf, t,
                                               anm[m + 1], pnm1);
#endif


                        /* Cm+1,m; Sm+1,m */
                        shcs->c[m][1] += SUM_R(MUL_R(pnm1, amm));
                        shcs->s[m][1] += SUM_R(MUL_R(pnm1, bmm));


                        /* Loop over degrees */
                        /* ------------------------------------------------- */
                        ds = 0;


                        /* Is "n + m" even?  Since we start the loop with "n
                         * = m + 2", then the parity of the first "m + 2 + m"
                         * is always even.  Then, it changes with every loop
                         * iteration. */
                        npm_even = 1;


                        for (unsigned long n = (m + 2); n <= nmax;
                             n++, npm_even = !npm_even)
                        {
                            /* Compute tesseral Legendre function */
#ifdef SIMD
                            PNM_TESSERAL_XNUM_SIMD(x, y, z, ix, iy, iz, ixy,
                                                   wlf, t, anm[n], bnm[n],
                                                   pnm2,
                                                   tmp1_r, tmp2_r,
                                                   mask1, mask2,
                                                   mask3, zero_r,
                                                   zero_ri, one_ri,
                                                   BIG_r, BIGI_r,
                                                   BIGS_r, BIGSI_r,
                                                   TESSERALS1, TESSERALS2, ds);
#else
                            PNM_TESSERAL_XNUM(x, y, z,
                                              ix, iy, iz,
                                              ixy, wlf, t,
                                              anm[n], bnm[n], pnm2, continue,
                                              ds);
#endif


                            /* Cm+2,m, Cm+3,m, ..., Cnmax,m and Sm+2,m, Sm+3,m,
                             * ..., Snmax,m */
                            nmm = n - m;
                            if (npm_even)
                            {
                                shcs->c[m][nmm] += SUM_R(MUL_R(pnm2, amp));
                                shcs->s[m][nmm] += SUM_R(MUL_R(pnm2, bmp));
                            }
                            else
                            {
                                shcs->c[m][nmm] += SUM_R(MUL_R(pnm2, amm));
                                shcs->s[m][nmm] += SUM_R(MUL_R(pnm2, bmm));
                            }


                        } /* End of the loop over harmonic degrees */
                        /* ------------------------------------------------- */


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
    REAL c2 = PREC(1.0) / (PREC(4.0) * PI) * (r0 / shcs->mu);
#if CHARM_OPENMP
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

