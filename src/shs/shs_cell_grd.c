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
#include "../leg/leg_func_anm_bnm.h"
#include "../leg/leg_func_dm.h"
#include "../leg/leg_func_gm_hm.h"
#include "../leg/leg_pol_en_fn.h"
#include "../leg/leg_func_r_ri.h"
#include "../leg/leg_func_prepare.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "shs_rpows.h"
#include "shs_check_symmi.h"
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


    /* Get the longitudinal step of the grid (will be necessary later for the
     * PSLR algorithm) */
    REAL dlon;
    if (cell_nlon > 1)
    {
        if (CHARM(misc_is_nearly_equal)(cell->lon[2] - cell->lon[0],
                                        cell->lon[3] - cell->lon[1],
                                        CHARM(glob_threshold2)) == 0)
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


    /* If true, FFT is applied along the latitude parallels */
    _Bool use_fft;


    /* Auxiliary constant to be used only in case the PSLR algorithm is
     * applied.  To suppress, a compiler warning, the value is initialized to
     * zero. */
    REAL lon0 = ADDP(0.0);


    /* Length of the lumped coefficients arrays in case FFT will be applied */
    size_t nlc = cell_nlon / 2 + 1;


    if (cell_nlon > 1)
        /* Let's check whether FFT can be employed.  Six conditions must be
         * satisfied to allow FFT:
         *
         * i) the number of grid longitudes must be large enough when compared
         * with "nmax",
         *
         * ii) the longitudinal step must be constant,
         *
         * iii) "cell->lon[0]" must be zero,
         *
         * iv) "cell->lon[cell_nlon 1] + dlon" must be equal to "2.0 * PI",
         *
         * v) there must be at least two cells in the longitudinal direction,
         * and
         *
         * vi) the maximum longitude of the "j"th cell must be equal to the
         * minimum longitude of the "j + 1"th cell.
         *
         * The second condition has already been checked, so let's do the
         * remaining checks.  Due to the previous check, the sixth condition
         * needs to be checked only for the first to longitude cells.*/
        if ((cell_nlon - 1) / 2 >= nmax &&
            CHARM(misc_is_nearly_equal)(cell->lon[0], ADDP(0.0),
                                        CHARM(glob_threshold)) &&
            CHARM(misc_is_nearly_equal)(cell->lon[2 * cell_nlon - 1],
                                        ADDP(2.0) * PI,
                                        CHARM(glob_threshold)) &&
            CHARM(misc_is_nearly_equal)(cell->lon[2], cell->lon[1],
                                        CHARM(glob_threshold)))
            /* Great, FFT can be applied for this grid! */
            use_fft = 1;
        else
            /* Oh no, FFT cannot be applied for this grid. */
        use_fft = 0;
    else
        /* There is only one cell in the longitudinal direction, so no FFT for
         * simplicity.  In case, the speed is not at all crucial. */
        use_fft = 0;


    if (!use_fft)
    {
        /* OK, so no FFT and the longitudinal step is constant */


        /* Get the origin of the longitude "lon" array (will be necessary later
         * for the PSLR algorithm).  The "lon" array is not actually stored in
         * memory, but is given as "(cell->lon[2 * j] + cell->lon[2 * j + 1])
         * / 2.0" for "j = 0, 1, ..., nlon - 1". */
        lon0 = (cell->lon[0] + cell->lon[1]) / ADDP(2.0);
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
        lc  = (FFTW(complex) *)FFTW(malloc)(nlc * sizeof(FFTW(complex)));
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


        int   *ips1        = NULL;
        int   *ips2        = NULL;
        REAL *ps1          = NULL;
        REAL *ps2          = NULL;
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


        ips1 =    (int *)calloc(nmax, sizeof(int));
        if (ips1 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        ips2 =    (int *)calloc(nmax, sizeof(int));
        if (ips2 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        ps1  = (REAL *)calloc(nmax, sizeof(REAL));
        if (ps1 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        ps2  = (REAL *)calloc(nmax, sizeof(REAL));
        if (ps2 == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        anm  = (REAL *)calloc(nmax + 1, sizeof(REAL));
        if (anm == NULL)
        {
            FAILURE_priv = 1;
            goto FAILURE_1_parallel;
        }


        bnm  = (REAL *)calloc(nmax + 1, sizeof(REAL));
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


            lc = (FFTW(complex) *)FFTW(malloc)(nlc * sizeof(FFTW(complex)));
            if (lc == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
            memset(lc, 0, nlc * sizeof(FFTW(complex)));
        }
        else
        {
            fi = (REAL *)calloc(cell_nlon, sizeof(REAL));
            if (fi == NULL)
            {
                FAILURE_priv = 1;
                goto FAILURE_1_parallel;
            }
        }


        rpows = (REAL *)calloc(nmax + 1, sizeof(REAL));
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


                lc2 = (FFTW(complex) *)FFTW(malloc)(nlc *
                                                    sizeof(FFTW(complex)));
                if (lc2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
                memset(lc2, 0, nlc * sizeof(FFTW(complex)));
            }
            else
            {
                fi2 = (REAL *)calloc(cell_nlon, sizeof(REAL));
                if (fi2 == NULL)
                {
                    FAILURE_priv = 1;
                    goto FAILURE_1_parallel;
                }
            }


            rpows2 = (REAL *)calloc(nmax + 1, sizeof(REAL));
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
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                               CHARM_ERR_MALLOC_FAILURE);
            }


            /* OK, and now all threads go to the "FAILURE_2_parallel" label
             * to deallocate all the memory that might be allocated before the
             * allocation failure. */
            goto FAILURE_2_parallel;
        }
        /* ................................................................. */


        REAL latmini, latmaxi;
        REAL t1, u1;
        REAL t2, u2;


        REAL a, b, a2, b2;
        a = b = a2 = b2 = ADDP(0.0);
        REAL c, cm, sm, mr;
        REAL imm0, imm1, imm2;
        REAL lontmp, clontmp, slontmp, dm0, dm1, dm2, dm02, dm12, dm22;
        dm02 = dm12 = dm22 = ADDP(0.0);
        REAL m2, m2sm2dl, cmdlon2;
        REAL dsigma;
        REAL mur_dsigma;

        _Bool symmi;
        size_t row_idx;


#if CHARM_PARALLEL
#pragma omp for
#endif
        for (size_t i = 0; i < nlatdo; i++)
        {
            /* Check whether the symmetry property of LFs needs to be exploited
             * */
            /* ------------------------------------------------------------- */
            symmi = CHARM(shs_check_symmi)(cell_type, nlatdo, symm, even, i);
            /* ------------------------------------------------------------- */


            /* Pre-compute the powers of "shcs->r / cell->r[i]" */
            /* ------------------------------------------------------------- */
            CHARM(shs_rpows)(shcs->r, cell->r[i], rpows, nmax);


            if (symmi)
                CHARM(shs_rpows)(shcs->r, cell->r[cell_nlat - i - 1],
                                 rpows2, nmax);
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


            /* ------------------------------------------------------------- */
            if (use_fft)
            {
                /* Important note.  For the "c2r" FFT transform, FFTW
                 * overwrites even the *input* array.  For some combinations of
                 * "nmax" and "cell_nlon", this is not a problem, while for
                 * other combinations, this causes incorrect results.  Here, we
                 * have to therefore reset "lc" and "lc2" to zeros.*/
                memset(lc, 0, nlc * sizeof(FFTW(complex)));


                if (symmi)
                    memset(lc2, 0, nlc * sizeof(FFTW(complex)));
            }
            else
            {
                /* The "fi" vector represents the synthesized quantity "f" for
                 * the "i"th latitude parallel. Therefore, it needs to be
                 * reinitialized to zero for each "ith" latitude. The same
                 * holds true for "fi2" in case of symmetric grids. */
                memset(fi, 0, cell_nlon * sizeof(REAL));


                if (symmi)
                    memset(fi2, 0, cell_nlon * sizeof(REAL));
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
                                       latmini, latmaxi,
                                       t1, t2,
                                       u1, u2,
                                       ps1, ps2,
                                       ips1, ips2,
                                       &imm0, &imm1, &imm2,
                                       en, fn, gm, hm, ri, rpows, rpows2,
                                       symmi,
                                       &a, &b, &a2, &b2);


                if (use_fft)
                {
                    c   = (m == 0) ? ADDP(1.0) : ADDP(0.5);


                    /* Due to the use of the mean values, some additional terms
                     * need to be taken into account when compared with the
                     * harmonic analysis based on point values (see, e.g.,
                     * Colombo, 1981) */
                    /* ----------------------------------------------------- */
                    if (m == 0)
                    {
                        sm = dlon;
                        cm = ADDP(0.0);
                    }
                    else
                    {
                        /* Useful substitution */
                        mr = (REAL)m;
                        cm = (COS(mr * dlon) - ADDP(1.0)) / mr;
                        sm = SIN(mr * dlon) / mr;
                    }


                    lc[m][0] =  (a * sm - b * cm) * c;
                    lc[m][1] = -(a * cm + b * sm) * c;

                    if (symmi)
                    {
                        lc2[m][0] =  (a2 * sm - b2 * cm) * c;
                        lc2[m][1] = -(a2 * cm + b2 * sm) * c;
                    }
                    /* ----------------------------------------------------- */
                }
                else
                {
                    /* Before applying the PSLR algorithm from Balmino et
                     * al. (2012), * the lumped coefficients need to be
                     * multiplied by additional terms, which stem from the
                     * integration in the longitudinal direction */
                    /* ----------------------------------------------------- */
                    if (m == 0)
                        m2sm2dl = dlon;
                    else
                    {
                        m2 = ADDP(2.0) / (REAL)m;
                        m2sm2dl = m2 * SIN(dlon / m2);
                    }


                    a *= m2sm2dl;
                    b *= m2sm2dl; /* Remember that "b = 0.0" for "m == 0.0", so
                                   * it can safely be multiplied by "m2sm2dl
                                   * * = dlon" */

                    if (symmi)
                    {
                        a2 *= m2sm2dl;
                        b2 *= m2sm2dl;  /* Remember that "b2 = 0.0" for "m ==
                                         * 0.0", so it can safely be multiplied
                                         * by "m2sm2dl = dlon" */
                    }
                    /* ----------------------------------------------------- */


                    /* The PSLR algorithm from Balmino et al. (2012) (see the
                     * REFERENCES section at the beginning of this function) */
                    /* ----------------------------------------------------- */

                    /* The first longitude cell */
                    /* ..................................................... */
                     lontmp = (REAL)m * lon0;
                    clontmp = COS(lontmp);
                    slontmp = SIN(lontmp);
                    dm0     = a * clontmp + b * slontmp;
                    fi[0]  += dm0;


                    if (symmi)
                    {
                        dm02    = a2 * clontmp + b2 * slontmp;
                        fi2[0] += dm02;
                    }


                    if (cell_nlon == 1)
                        continue;
                    /* ..................................................... */


                    /* The second longitude cell */
                    /* ..................................................... */
                     lontmp = (REAL)m * (lon0 + dlon);
                    clontmp = COS(lontmp);
                    slontmp = SIN(lontmp);
                    dm1     = a * clontmp + b * slontmp;
                    fi[1]  += dm1;


                    if (symmi)
                    {
                        dm12    = a2 * clontmp + b2 * slontmp;
                        fi2[1] += dm12;
                    }


                    if (cell_nlon == 2)
                        continue;
                    /* ..................................................... */


                    /* The third and all the remaining longitude cells */
                    /* ..................................................... */
                    cmdlon2 = ADDP(2.0) * COS((REAL)m * dlon);
                    for (size_t j = 2; j < cell_nlon; j++)
                    {
                        dm2    = cmdlon2 * dm1 - dm0;
                        fi[j] += dm2;
                        dm0    = dm1;
                        dm1    = dm2;
                    }


                    if (symmi)
                    {
                        for (size_t j = 2; j < cell_nlon; j++)
                        {
                            dm22    = cmdlon2 * dm12 - dm02;
                            fi2[j] += dm22;
                            dm02    = dm12;
                            dm12    = dm22;
                        }
                    }
                    /* ..................................................... */
                    /* ----------------------------------------------------- */
                    /* End of the PSLR algorithm */
                }

            } /* End of the loop over harmonic orders */
            /* ------------------------------------------------------------- */


            /* Cell area on the unit sphere and some useful constants */
            dsigma = (SIN(latmaxi) - SIN(latmini)) * dlon;
            mur_dsigma = mur / dsigma;
            row_idx = i * cell_nlon;


            if (use_fft)
            {
                /* Fourier transform along the latitude parallels */
                /* --------------------------------------------------------- */
                FFTW(execute_dft_c2r)(plan, lc, ftmp);
                for (size_t j = 0; j < cell_nlon; j++)
                    f[row_idx + j] = mur_dsigma * ftmp[j];


                if (symmi)
                {
                    row_idx = (cell_nlat - i - 1) * cell_nlon;
                    FFTW(execute_dft_c2r)(plan, lc2, ftmp2);
                    for (size_t j = 0; j < cell_nlon; j++)
                        f[row_idx + j] = mur_dsigma * ftmp2[j];
                }
                /* --------------------------------------------------------- */
            }
            else
            {
                /* Final synthesis */
                /* --------------------------------------------------------- */
                for (size_t j = 0; j < cell_nlon; j++)
                    f[row_idx + j] = mur_dsigma * fi[j];


                if (symmi)
                {
                    row_idx = (cell_nlat - i - 1) * cell_nlon;
                    for (size_t j = 0; j < cell_nlon; j++)
                        f[row_idx + j] = mur_dsigma * fi2[j];
                }
                /* --------------------------------------------------------- */
            }


        } /* End of the loop over latitude parallels */
        /* ----------------------------------------------------------------- */


FAILURE_2_parallel:
        free(ips1); free(ps1);
        free(ips2); free(ps2);
        free(anm); free(bnm);
        free(fi); free(fi2);
        FFTW(free)(lc); FFTW(free)(lc2);
        FFTW(free)(ftmp); FFTW(free)(ftmp2);
        free(rpows); free(rpows2);


    } /* End of "#pragma omp parallel" */
    /* --------------------------------------------------------------------- */






    /* Freeing up the heap memory */
    /* --------------------------------------------------------------------- */
FAILURE:
    if (FAILURE_glob != 0)
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
