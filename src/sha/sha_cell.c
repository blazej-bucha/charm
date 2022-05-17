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
#include "../leg/leg_func_gm_hm.h"
#include "../leg/leg_pol_en_fn.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../leg/leg_func_xnum.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#if CHARM_PARALLEL
#   include <omp.h>
#endif
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
    /* ..................................................................... */


    if (CHARM(misc_is_nearly_equal)(cell->lon[0], ADDP(0.0),
                                    CHARM(glob_threshold)) != 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lon[0]\" has to be equal to \"0.0\" "
                       "within the \"threshold\".");
        return;
    }


    if (CHARM(misc_is_nearly_equal)(cell->lon[2 * cell->nlon - 1],
                                    ADDP(2.0) * PI,
                                    CHARM(glob_threshold)) != 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lon[2 * cell->nlon - 1]\" has to be "
                       "equal to \"2.0 * PI\" within the "
                       "\"threshold\".");
        return;
    }


    /* Check whether "cell->lon" is a linearly increasing array of cells. */
    int err_tmp = CHARM(misc_arr_chck_lin_incr)(cell->lon, 2 * cell_nlon,
                                                0, 2, CHARM(glob_threshold2),
                                                err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    if (err_tmp != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lon\" is not a linearly increasing array "
                       "of cells within the \"threshold2\".");
        return;
    }


    err_tmp = CHARM(misc_arr_chck_lin_incr)(cell->lon, 2 * cell_nlon,
                                            1, 2, CHARM(glob_threshold2), err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    if (err_tmp != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cell->lon\" is not a linearly increasing array "
                       "of cells within the \"threshold2\".");
        return;
    }


    /* At this point, we know that the steps between cells in "cell->lon" are
     * constant.  Here, we check if the steps between the minimum longitudes
     * and the maximum longitudes are equal and, if true, store this value in
     * a separate variable (will be necessary later for the PSLR algorithm) */
    REAL dlon;
    if (cell_nlon > 1)
    {
        if (CHARM(misc_is_nearly_equal)(cell->lon[2] - cell->lon[0],
                                        cell->lon[3] - cell->lon[1],
                                        CHARM(glob_threshold)) == 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The difference \"cell->lon[2] - cell->lon[0]\" "
                           "has to be equal to "
                           "\"cell->lon[3] - cell->lon[1]\".");


            return;
        }
        dlon = cell->lon[2] - cell->lon[0];
    }
    else
        dlon = ADDP(0.0);


    /* Now check whether the grid cells are nicely aligned in the longitudinal
     * direction.  Due to the previous error checks, this can be done
     * only with the first grid cell in the longitudinal direction. */
    if (CHARM(misc_is_nearly_equal)(cell->lon[2], cell->lon[1],
                                    CHARM(glob_threshold)) == 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The grid cells are not nicely aligned in the "
                       "longitudinal direction.  It must hold that "
                       "\"cell->lon[2 * j] == cell->lon[2 * j - 1]\" for all "
                       "\"j = 1, 2, ..., cell->nlon - 1\".");
        return;
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
    REAL *ftmp_in        = NULL;
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






    /* Initialize all elements of the output coefficients to zero */
    /* --------------------------------------------------------------------- */
    unsigned long nm_count = ((nmax + 2) * (nmax + 1)) / 2;
    memset(shcs->c[0], 0, nm_count * sizeof(REAL));
    memset(shcs->s[0], 0, nm_count * sizeof(REAL));
    /* --------------------------------------------------------------------- */






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
    ftmp_out = (FFTW(complex) *)FFTW(malloc)(sizeof(FFTW(complex))
                                             * (cell_nlon / 2 + 1));
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
    size_t imax;
    if (symm)
        imax = nlatdo + even - 1;
    else
        /* The grid is not symmetric, so we do not want to apply the symmetry
         * property of Legendre functions. */
        imax = 0;


    {
        REAL latmini, latmaxi;
        REAL t1, u1, x1, y1, z1;
        REAL t2, u2, x2, y2, z2;
        int ix1, iy1, iz1, ixy1;
        int ix2, iy2, iz2, ixy2;
        REAL pnm0_latmini, pnm1_latmini, pnm2_latmini;
        REAL pnm0_latmaxi, pnm1_latmaxi, pnm2_latmaxi;
        REAL w;
        REAL amp, amm, bmp, bmm;
        REAL cm, sm;
        REAL in0, inm0, inm1, inm2;
        REAL mr;
        _Bool symmi;
        _Bool npm_even; /* True if "n + m" is even */


        /* ................................................................. */
        int   *ips1 = NULL;
        int   *ips2 = NULL;
        REAL *ps1 = NULL;
        REAL *ps2 = NULL;
        REAL *a   = NULL;
        REAL *b   = NULL;
        REAL *a2  = NULL;
        REAL *b2  = NULL;
#if !(CHARM_PARALLEL)
        REAL *anm = NULL;
        REAL *bnm = NULL;
#endif
        REAL *imm = NULL;
        REAL *ftmp_in = NULL;
        FFTW(complex) *ftmp_out = NULL;


        ips1 = (int *)calloc(nmax + 1, sizeof(int));
        if (ips1 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ips2 = (int *)calloc(nmax + 1, sizeof(int));
        if (ips2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ps1 = (REAL *)calloc(nmax, sizeof(REAL));
        if (ps1 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        ps2 = (REAL *)calloc(nmax, sizeof(REAL));
        if (ps2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        a   = (REAL *)calloc(cell_nlon / 2 + 1, sizeof(REAL));
        if (a == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b   = (REAL *)calloc(cell_nlon / 2 + 1, sizeof(REAL));
        if (b == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        a2  = (REAL *)calloc(cell_nlon / 2 + 1, sizeof(REAL));
        if (a2 == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        b2  = (REAL *)calloc(cell_nlon / 2 + 1, sizeof(REAL));
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
        imm = (REAL *)calloc((nmax + 1), sizeof(REAL));
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
                                                 (cell_nlon / 2 + 1));
        if (ftmp_out == NULL)
        {
            FAILURE_glob = 1;
            goto FAILURE_1;
        }
        /* ................................................................. */



        for (size_t i = 0; i < nlatdo; i++)
        {

            /* Lumped coefficients for the southern hemisphere (including the
             * equator) */
            /* ------------------------------------------------------------- */
            memcpy(ftmp_in, f + i * cell_nlon, cell_nlon * sizeof(REAL));
            FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


            for (size_t j = 0; j < (cell_nlon / 2 + 1); j++)
            {
                a[j] =  ftmp_out[j][0];
                b[j] = -ftmp_out[j][1];
            }
            /* ------------------------------------------------------------- */


            /* Lumped coefficients for the northern hemisphere */
            /* ------------------------------------------------------------- */
            symmi = i < imax;
            if (symmi)
            {
                memcpy(ftmp_in, f + (cell_nlat - i - 1) * cell_nlon,
                       cell_nlon * sizeof(REAL));
                FFTW(execute_dft_r2c)(plan, ftmp_in, ftmp_out);


                for (size_t j = 0; j < (cell_nlon / 2 + 1); j++)
                {
                    a2[j] =  ftmp_out[j][0];
                    b2[j] = -ftmp_out[j][1];
                }
            }
            /* ------------------------------------------------------------- */


            /* Prepare arrays for sectorial Legendre functions */
            /* ------------------------------------------------------------- */
            latmini = cell->lat[2 * i];
            latmaxi = cell->lat[2 * i + 1];


            t1    = SIN(latmini);
            u1    = COS(latmini);
            CHARM(leg_func_prepare)(u1, ps1, ips1, dm, nmax);


            t2    = SIN(latmaxi);
            u2    = COS(latmaxi);
            CHARM(leg_func_prepare)(u2, ps2, ips2, dm, nmax);
            /* ------------------------------------------------------------- */


            /* Pre-compute the sectorial "imm" integrals.  This is necessary
             * for the "defined(CHARM_PARALLEL)" parallelization strategy, but
             * can be used (and in fact it is) also with the other
             * parallelization strategies. */
            /* ------------------------------------------------------------- */
            for (unsigned long m = 1; m <= nmax; m++)
            {
                /* Pmm for latmini */
                PNM_SECTORIAL_XNUM(x1, ix1, ps1[m - 1], ips1[m - 1],
                                   pnm0_latmini);


                /* Pmm for latmaxi */
                PNM_SECTORIAL_XNUM(x2, ix2, ps2[m - 1], ips2[m - 1],
                                   pnm0_latmaxi);


                /* Imm */
                /* ..................................................... */
                if (m == 1) /* I11 */
                {
                    imm[m] = SQRT(ADDP(3.0)) / ADDP(2.0) *
                             ((t2 * u2 - (PI_2 - latmaxi)) -
                              (t1 * u1 - (PI_2 - latmini)));
                }
                else if (m == 2) /* I22 */
                {
                    imm[m] = SQRT(ADDP(15.0)) / ADDP(6.0) *
                             (t2 * (ADDP(3.0) - t2 * t2) -
                              t1 * (ADDP(3.0) - t1 * t1));
                }
                else /* I33, I44, ..., Inmax,nmax */
                {
                    imm[m] = gm[m] * imm[m - 2] + ADDP(1.0) / (REAL)(m + 1) *
                             (t2 * pnm0_latmaxi - t1 * pnm0_latmini);
                }
                /* ..................................................... */
            }
            /* ------------------------------------------------------------- */






            /* Loop over harmonic orders */
            /* ------------------------------------------------------------- */
#if CHARM_PARALLEL
#pragma omp parallel default(none) \
shared(nmax, t1, t2, u1, u2, latmini, latmaxi, i, symmi) \
shared(shcs, en, fn, gm, hm, imm, ps1, ps2, ips1, ips2, r, ri) \
shared(a, b, a2, b2, dlon) \
shared(FAILURE_glob, err) \
private(cm, sm, x1, x2, ix1, ix2, y1, y2, iy1, iy2, w, mr) \
private(z1, z2, iz1, iz2, ixy1, ixy2) \
private(amp, amm, bmp, bmm) \
private(pnm0_latmini, pnm0_latmaxi, pnm1_latmini, pnm1_latmaxi) \
private(pnm2_latmini, pnm2_latmaxi, in0, inm0, inm1, inm2, npm_even)
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
                {
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EMEM, CHARM_ERR_MALLOC_FAILURE);
                }


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
                if (m == 0)
                {
                    if (symmi)
                    {
                        amp = (a[m] + a2[m]) * dlon;
                        amm = (a[m] - a2[m]) * dlon;
                    }
                    else
                    {
                        amp = a[m] * dlon;
                        amm = amp;
                    }
                }
                else
                {
                    /* Useful substitution */
                    mr = (REAL)m;
                    cm = (COS(mr * dlon) - ADDP(1.0)) / mr;
                    sm = SIN(mr * dlon) / mr;


                    if (symmi)
                    {
                        amp =  (a[m] + a2[m]) * sm + (b[m] + b2[m]) * cm;
                        amm =  (a[m] - a2[m]) * sm + (b[m] - b2[m]) * cm;
                        bmp = -(a[m] + a2[m]) * cm + (b[m] + b2[m]) * sm;
                        bmm = -(a[m] - a2[m]) * cm + (b[m] - b2[m]) * sm;
                    }
                    else
                    {
                        amp = a[m] * sm + b[m] * cm;
                        amm = amp;
                        bmp = -a[m] * cm + b[m] * sm;
                        bmm = bmp;
                    }
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
                    pnm0_latmini = ADDP(1.0);
                    pnm0_latmaxi = ADDP(1.0);
                    /* ..................................................... */


                    /* P10 */
                    /* ..................................................... */
                    pnm1_latmini = t1;
                    pnm1_latmaxi = t2;
                    /* ..................................................... */


                    /* I00 */
                    /* ..................................................... */
                    in0 = t2 - t1;
                    /* ..................................................... */


                    /* Spherical harmonic coefficients */
                    /* ..................................................... */
                    shcs->c[0][0] += in0 * amp;
                    /* ..................................................... */


                    if ((nmax + 1) >= 2)
                    {
                        for (unsigned long n = 1; n <= nmax; n++)
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
                            pnm2_latmini = en[n + 1] * t1 * pnm1_latmini -
                                           fn[n + 1] * pnm0_latmini;
                            pnm2_latmaxi = en[n + 1] * t2 * pnm1_latmaxi -
                                           fn[n + 1] * pnm0_latmaxi;
                            /* ............................................. */


                            /* I10, I20, ..., Inmax,0
                             * Computed from Pn+1,0 and Pn-1,0 */
                            /* ............................................. */
                            in0 = ri[2 * n + 1] *
                                  (pnm2_latmaxi - pnm0_latmaxi -
                                   pnm2_latmini + pnm0_latmini);
                            /* ............................................. */


                            /* Spherical harmonic coefficients */
                            /* ............................................. */
                            /* C10, C20, ..., Cnmax,0 */
                            if ((n % 2) == 0)
                                shcs->c[0][n] += in0 * amp;
                            else
                                shcs->c[0][n] += in0 * amm;
                            /* ............................................. */


                            pnm0_latmini = pnm1_latmini;
                            pnm1_latmini = pnm2_latmini;


                            pnm0_latmaxi = pnm1_latmaxi;
                            pnm1_latmaxi = pnm2_latmaxi;

                        }
                    }
                    /* ----------------------------------------------------- */







                }
                else /* Non-zonal harmonics */
                {






                    /* Sectorial Legendre functions and their integrals */
                    /* ----------------------------------------------------- */

                    /* Pmm for latmini */
                    PNM_SECTORIAL_XNUM(x1, ix1, ps1[m - 1], ips1[m - 1],
                                       pnm0_latmini);


                    /* Pmm for latmaxi */
                    PNM_SECTORIAL_XNUM(x2, ix2, ps2[m - 1], ips2[m - 1],
                                       pnm0_latmaxi);


                    /* Imm */
                    /* ..................................................... */
                    inm0 = imm[m];
                    /* ..................................................... */


                    /* Spherical harmonic coefficients */
                    /* ..................................................... */
                    shcs->c[m][0] += inm0 * amp;
                    shcs->s[m][0] += inm0 * bmp;
                    /* ..................................................... */
                    /* ----------------------------------------------------- */






                    /* Tesseral harmonics */
                    /* ----------------------------------------------------- */
                    if (m < nmax)
                    {

                        /* Pm+1,m for latmini */
                        PNM_SEMISECTORIAL_XNUM(x1, y1, ix1, iy1, w, t1,
                                               anm[m + 1], pnm1_latmini);


                        /* Pm+1,m for latmaxi */
                        PNM_SEMISECTORIAL_XNUM(x2, y2, ix2, iy2, w, t2,
                                               anm[m + 1], pnm1_latmaxi);


                        /* Im+1,m */
                        /* ................................................. */
                        /* This is not a typo, "pnm0" are indeed required here
                         * */
                        inm1 = -anm[m + 1] / (REAL)(m + 2) *
                                (u2 * u2 * pnm0_latmaxi -
                                 u1 * u1 * pnm0_latmini);
                        /* ................................................. */


                        /* Spherical harmonic coefficients */
                        /* ................................................. */
                        shcs->c[m][1] += inm1 * amm;
                        shcs->s[m][1] += inm1 * bmm;
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

                            /* Pm+2,m, Pm+3,m, ..., Pnmax,m for latmini */
                            PNM_TESSERAL_XNUM(x1, y1, z1,
                                              ix1, iy1, iz1,
                                              ixy1, w, t1,
                                              anm[n], bnm[n],
                                              pnm2_latmini,
                                              pnm2_latmini = ADDP(0.0));


                            /* Pm+2,m, Pm+3,m, ..., Pnmax,m for latmaxi */
                            PNM_TESSERAL_XNUM(x2, y2, z2,
                                              ix2, iy2, iz2,
                                              ixy2, w, t2,
                                              anm[n], bnm[n],
                                              pnm2_latmaxi,
                                              pnm2_latmaxi = ADDP(0.0));


                            /* Im+2,m, Im+3,m, ..., Inmax,m */
                            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                            /* This is not a typo, "pnm1" are indeed required
                             * here */
                            inm2 = hm[n] * bnm[n] * inm0 -
                                   anm[n] / (REAL)(n + 1) *
                                   (u2 * u2 * pnm1_latmaxi -
                                    u1 * u1 * pnm1_latmini);

                            inm0 = inm1;
                            inm1 = inm2;
                            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


                            /* Spherical harmonic coefficients */
                            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                            /* Cm+2,m, Cm+3,m, ..., Cnmax,m and Sm+2,m, Sm+3,m,
                             * * ..., Snmax,m */
                            if (npm_even)
                            {
                                shcs->c[m][n - m] += inm2 * amp;
                                shcs->s[m][n - m] += inm2 * bmp;
                            }
                            else
                            {
                                shcs->c[m][n - m] += inm2 * amm;
                                shcs->s[m][n - m] += inm2 * bmm;
                            }
                            /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


                            pnm1_latmini = pnm2_latmini;
                            pnm1_latmaxi = pnm2_latmaxi;

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
        if (FAILURE_glob != 0)
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);


        FFTW(free)(ftmp_in);
        FFTW(free)(ftmp_out);
        free(ips1); free(ps1);
        free(ips2); free(ps2);
#if !(CHARM_PARALLEL)
        free(anm);  free(bnm);
#endif
        free(a);    free(b);
        free(a2);   free(b2);
        free(imm);


    }
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if (FAILURE_glob != 0)
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
