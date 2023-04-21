/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/prec.h"
#include "generate_point.h"
#include "parameters.h"
#include "validate.h"
/* ------------------------------------------------------------------------- */






long int check_shs_point(void)
{
    /* Read reference potential coefficients */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_pot = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_pot == NULL)
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


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs_pot, err);
    CHARM(err_handler)(err, 1);


    /* Modify coefficients of degrees "0" and "1" to allow for an accurate
     * validation in all precisions. */
    shcs_pot->c[0][0 - 0] = (REAL)C00;
    shcs_pot->c[0][1 - 0] = (REAL)C10;
    shcs_pot->c[1][1 - 1] = (REAL)C11;
    shcs_pot->s[1][1 - 1] = (REAL)S11;
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
                REAL rref = shcs_pot->r + (REAL)(DELTAR) * (REAL)deltar;


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


                f = (REAL *)malloc(grd_pnt->nlat * grd_pnt->nlon *
                                   sizeof(REAL));
                if (f == NULL)
                {
                    fprintf(stderr, "malloc failure.\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(shs_point)(grd_pnt, shcs_pot, nmax, f, err);
                CHARM(err_handler)(err, 1);


                e += validate(file, f, grd_pnt->nlat * grd_pnt->nlon,
                              PREC(10.0) * CHARM(glob_threshold));


                CHARM(crd_point_free)(grd_pnt);
                free(f);
            }
        }
    }
    }
    /* ..................................................................... */






    /* Custom grid of points and cells */
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
                        REAL r = shcs_pot->r + (REAL)(DELTAR) * (REAL)deltar;


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


                        f = (REAL *)malloc(grd_pnt->nlat * grd_pnt->nlon *
                                           sizeof(REAL));
                        if (f == NULL)
                        {
                            fprintf(stderr, "malloc failure.\n");
                            exit(CHARM_FAILURE);
                        }


                        CHARM(shs_point)(grd_pnt, shcs_pot, nmax, f, err);
                        CHARM(err_handler)(err, 1);
                        e += validate(file, f,
                                      grd_pnt->nlat *
                                      grd_pnt->nlon,
                                      PREC(100.0) *
                                      CHARM(glob_threshold));


                        CHARM(crd_point_free)(grd_pnt);


                        free(f);
                    }
                }
            }
        }
    }
    }
    /* ..................................................................... */






    /* Scattered points/cells */
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
                REAL r = shcs_pot->r + (REAL)(DELTAR) * (REAL)deltar;


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


                f = (REAL *)malloc(sctr_pnt->nlat * sizeof(REAL));
                if (f == NULL)
                {
                    fprintf(stderr, "malloc failure.\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(shs_point)(sctr_pnt, shcs_pot, nmax, f, err);
                CHARM(err_handler)(err, 1);
                e += validate(file, f, sctr_pnt->nlat,
                              PREC(10.0) * CHARM(glob_threshold));


                CHARM(crd_point_free)(sctr_pnt);


                free(f);
            }
        }
    }
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);
    CHARM(shc_free)(shcs_pot);


    return e;
    /* --------------------------------------------------------------------- */
}
