/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/prec.h"
#include "generate_point.h"
#include "parameters.h"
#ifdef GENREF
#   include "write.h"
#else
#   include "validate.h"
#endif
/* ------------------------------------------------------------------------- */






long int check_shs_point(void)
{
    /* Read reference potential coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0), PREC(1.0));
    if (shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Error structure */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs, err);
    CHARM(err_handler)(err, 1);


    /* Modify coefficients of degrees "0" and "1" to allow for an accurate
     * validation in all precisions. */
    shcs->c[0][0 - 0] = (REAL)C00;
    shcs->c[0][1 - 0] = (REAL)C10;
    shcs->c[1][1 - 1] = (REAL)C11;
    shcs->s[1][1 - 1] = (REAL)S11;
    /* --------------------------------------------------------------------- */






    /* GL, DH1 and DH2 point grids */
    /* ..................................................................... */
    long int e            = 0;
    REAL *f               = NULL;
    CHARM(point) *grd_pnt = NULL;
    char file[NSTR_LONG];
    char grd_str[NSTR_SHORT];


    {
    int grd_types[3] = {CHARM_CRD_POINT_GRID_GL,
                        CHARM_CRD_POINT_GRID_DH1,
                        CHARM_CRD_POINT_GRID_DH2};


    for (int g = 0; g < 3; g++)
    {
        int grd_type = grd_types[g];


        for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
        {
            for (int deltar = 0; deltar < NDELTAR; deltar++)
            {
                REAL rref = shcs->r + (REAL)(DELTAR) * (REAL)deltar;


                if (grd_type == CHARM_CRD_POINT_GRID_GL)
                    grd_pnt = CHARM(crd_point_gl)(nmax, rref);
                else if (grd_type == CHARM_CRD_POINT_GRID_DH1)
                    grd_pnt = CHARM(crd_point_dh1)(nmax, rref);
                else if (grd_type == CHARM_CRD_POINT_GRID_DH2)
                    grd_pnt = CHARM(crd_point_dh2)(nmax, rref);
                if (grd_pnt == NULL)
                {
                    fprintf(stderr, "Failed to initialize a \"crd\" "
                            "structure\n");
                    exit(CHARM_FAILURE);
                }


                /* Generate output file name */
                if (g == 0)
                    strcpy(grd_str, "gl");
                else if (g == 1)
                    strcpy(grd_str, "dh1");
                else if (g == 2)
                    strcpy(grd_str, "dh2");
                sprintf(file, "%s/shs_p_nx%lu_dr%d_fft1_%s%s",
                        FOLDER, nmax, deltar, grd_str, FTYPE);


                f = (REAL *)malloc(grd_pnt->npoint * sizeof(REAL));
                if (f == NULL)
                {
                    fprintf(stderr, "malloc failure.\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(shs_point)(grd_pnt, shcs, nmax, f, err);
                CHARM(err_handler)(err, 1);


#ifdef GENREF
                e += write(file, f, grd_pnt->npoint);
#else
                e += validate(file, f, grd_pnt->npoint,
                              PREC(10.0) * CHARM(glob_threshold));
#endif


                CHARM(crd_point_free)(grd_pnt);
                free(f);
            }
        }
    }
    }
    /* ..................................................................... */






    /* Custom point grids */
    /* ..................................................................... */
    {
    size_t nlat[NCUSTOM_GRD] = {1, 1, 3, 10};
    size_t nlon[NCUSTOM_GRD] = {1, 2, 8, 22};


    for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
    {
        for (size_t i = 0; i < NCUSTOM_GRD; i++)
        {
            for (int fft = 0; fft < 2; fft++)
            {
                if (fft == 1)
                {
                    if (((nlon[i] - 1) / 2 < nmax) || (nlon[i] < 2))
                        continue;
                }


                for (int s = 0; s < 2; s++) /* "s = 0" for non-symm grids
                                             */
                {
                    if ((nlat[i] == 1) && (s == 1))
                        /* If there is only a single latitude, the grid is
                         * automatically non-symmetric */
                        continue;


                    for (int deltar = 0; deltar < NDELTAR; deltar++)
                    {
                        REAL r = shcs->r + (REAL)(DELTAR) * (REAL)deltar;


                        grd_pnt = CHARM(crd_point_calloc)(CHARM_CRD_POINT_GRID,
                                                          nlat[i], nlon[i]);
                        if (grd_pnt == NULL)
                        {
                            fprintf(stderr, "Failed to initialize a "
                                            "\"crd\" structure\n");
                            exit(CHARM_FAILURE);
                        }


                        if (fft == 0)
                            CHARM(generate_point)(grd_pnt, r, PI, PI);
                        else
                            CHARM(generate_point)(grd_pnt, r, PI,
                                                  PREC(2.0) * PI);


                        /* To get a non-symmetric grid, we simply add
                         * some more or less random number to the first
                         * latitude */
                        REAL break_symm = PREC(0.0);
                        if (s == 0)
                            break_symm = (REAL)(BREAK_SYMM);
                        grd_pnt->lat[0] -= break_symm;


                        /* Generate output file name */
                        sprintf(file, "%s/shs_%s_nx%lu_n%zu_dr%d_fft%d"
                                          "_s%d%s",
                                FOLDER, "p", nmax, i, deltar, fft,
                                (s == 0) ? 0 : 1, FTYPE);


                        f = (REAL *)malloc(grd_pnt->npoint * sizeof(REAL));
                        if (f == NULL)
                        {
                            fprintf(stderr, "malloc failure.\n");
                            exit(CHARM_FAILURE);
                        }


                        CHARM(shs_point)(grd_pnt, shcs, nmax, f, err);
                        CHARM(err_handler)(err, 1);
#ifdef GENREF
                        e += write(file, f, grd_pnt->npoint);
#else
                        e += validate(file, f, grd_pnt->npoint,
                                      PREC(100.0) * CHARM(glob_threshold));
#endif


                        CHARM(crd_point_free)(grd_pnt);


                        free(f);
                    }
                }
            }
        }
    }
    }
    /* ..................................................................... */






    /* Scattered points */
    /* ..................................................................... */
    {
    size_t nlat[NSCTR] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    size_t nlon[NSCTR] = {1, 2, 3, 4, 5, 6, 7, 8, 9};


    CHARM(point) *sctr_pnt = NULL;
    for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
    {
        for (size_t i = 0; i < NSCTR; i++)
        {
            for (int deltar = 0; deltar < NDELTAR; deltar++)
            {
                REAL r = shcs->r + (REAL)(DELTAR) * (REAL)deltar;


                sctr_pnt = CHARM(crd_point_malloc)(CHARM_CRD_POINT_SCATTERED,
                                                   nlat[i], nlon[i]);
                if (sctr_pnt == NULL)
                {
                    fprintf(stderr, "Failed to initialize a \"crd\" "
                            "structure\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(generate_point)(sctr_pnt, r, PI, PREC(2.0) * PI);


                /* Generate output file name */
                sprintf(file, "%s/shs_%s_nx%lu_n%zu_dr%d_sctr%s",
                        FOLDER, "p", nmax, i, deltar, FTYPE);


                f = (REAL *)malloc(sctr_pnt->npoint * sizeof(REAL));
                if (f == NULL)
                {
                    fprintf(stderr, "malloc failure.\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(shs_point)(sctr_pnt, shcs, nmax, f, err);
                CHARM(err_handler)(err, 1);
#ifdef GENREF
                e += write(file, f, sctr_pnt->npoint);
#else
                e += validate(file, f, sctr_pnt->npoint,
                              PREC(10.0) * CHARM(glob_threshold));
#endif


                CHARM(crd_point_free)(sctr_pnt);


                free(f);
            }
        }
    }
    }


    CHARM(shc_free)(shcs);
    /* --------------------------------------------------------------------- */






    /* Dynamical switching and loop unrolling at custom point grids.  Similarly
     * as in the tests of harmonic analysis, we use some fake coefficients.
     * See "check_sha_point.c" for further details. */
    /* ..................................................................... */
    {
    size_t nlat[NCUSTOM_GRD] = {10, 11, 12, 13};
    size_t nlon[NCUSTOM_GRD] = {22, 24, 26, 28};


    for (unsigned long nmax = NMAX_DS_MIN; nmax <= NMAX_DS_MAX; nmax++)
    {
        shcs = CHARM(shc_calloc)(nmax, PREC(1.0), PREC(1.0));
        if (shcs == NULL)
        {
            fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
            exit(CHARM_FAILURE);
        }
        for (size_t m = 0; m <= nmax; m++)
        {
            for (size_t n = m; n <= nmax; n++)
            {
                shcs->c[m][n - m] = PREC(1.0) / (REAL)(n + 1);
                if (m > 0)
                    shcs->s[m][n - m] = PREC(1.0) / (REAL)(n + 1);
            }
        }


        for (size_t i = 0; i < NCUSTOM_GRD; i++)
        {
            for (int fft = 0; fft < 2; fft++)
            {
                if (fft == 1)
                {
                    if (((nlon[i] - 1) / 2 < nmax) || (nlon[i] < 2))
                        continue;
                }


                for (int s = 0; s < 2; s++) /* "s = 0" for non-symm grids */
                {
                    if ((nlat[i] == 1) && (s == 1))
                        /* If there is only a single latitude, the grid is
                         * automatically non-symmetric */
                        continue;


                    for (int deltar = 0; deltar < NDELTAR; deltar++)
                    {
                        REAL r = shcs->r + (REAL)(DELTAR) * (REAL)deltar;


                        grd_pnt = CHARM(crd_point_gl)(nmax, r);
                        if (grd_pnt == NULL)
                        {
                            fprintf(stderr, "Failed to initialize a "
                                            "\"crd\" structure\n");
                            exit(CHARM_FAILURE);
                        }
                        /* Given that we are testing features that are related
                         * to latitudes only, we can do a nasty thing and
                         * modify the number of longitudes in "grd_pnt" to "1",
                         * so that the synthesis will be done only for the
                         * first meridian in "grd_pnt".  This makes the tests
                         * faster and avoids having too large reference data
                         * files in "../data/tests". */
                        grd_pnt->nlon = 1;


                        /* To get a non-symmetric grid, we simply add
                         * some more or less random number to the first
                         * latitude */
                        REAL break_symm = PREC(0.0);
                        if (s == 0)
                            break_symm = (REAL)(BREAK_SYMM);
                        grd_pnt->lat[0] -= break_symm;


                        /* Generate output file name */
                        sprintf(file, "%s/shs_%s_nx%lu_n%zu_dr%d_fft%d"
                                          "_s%d_dsun%s",
                                FOLDER, "p", nmax, i, deltar, fft,
                                (s == 0) ? 0 : 1, FTYPE);


                        f = (REAL *)malloc(grd_pnt->npoint * sizeof(REAL));
                        if (f == NULL)
                        {
                            fprintf(stderr, "malloc failure.\n");
                            exit(CHARM_FAILURE);
                        }


                        CHARM(shs_point)(grd_pnt, shcs, nmax, f, err);
                        CHARM(err_handler)(err, 1);
#ifdef GENREF
                        /* Do not use "grd_pnt->npoint" here, because
                         * "grd_pnt->nlon" has been modified. */
                        e += write(file, f, grd_pnt->nlat * grd_pnt->nlon);
#else
                        /* Do not use "grd_pnt->npoint" here, because
                         * "grd_pnt->nlon" has been modified. */
                        e += validate(file, f, grd_pnt->nlat * grd_pnt->nlon,
                                      CHARM(glob_threshold2));
#endif


                        CHARM(crd_point_free)(grd_pnt);


                        free(f);
                    }
                }
            }
        }
    }
    }
    /* ..................................................................... */






    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);
    CHARM(shc_free)(shcs);


    return e;
    /* --------------------------------------------------------------------- */
}
