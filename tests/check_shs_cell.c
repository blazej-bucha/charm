/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "generate_cell.h"
#include "parameters.h"
#ifdef GENREF
#   include "array2file.h"
#else
#   include "validate.h"
#endif
#include "modify_low_degree_coefficients.h"
#include "check_shs_cell.h"
/* ------------------------------------------------------------------------- */






long int check_shs_cell(void)
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
    modify_low_degree_coefficients(shcs_pot);
    /* --------------------------------------------------------------------- */






    /* Custom cell grids */
    /* ..................................................................... */
    long int e            = 0;
    REAL *f               = NULL;
    CHARM(cell) *grd_cell = NULL;
    char file[NSTR_LONG];


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
                        REAL r = shcs_pot->r + (REAL)(DELTAR) *
                                 (REAL)deltar;


                        grd_cell = CHARM(crd_cell_calloc)(CHARM_CRD_CELL_GRID,
                                                          nlat[i], nlon[i]);
                        if (grd_cell == NULL)
                        {
                            fprintf(stderr, "Failed to initialize a "
                                            "\"crd\" structure\n");
                            exit(CHARM_FAILURE);
                        }


                        if (fft == 0)
                            CHARM(generate_cell)(grd_cell, r, PI, PI);
                        else
                            CHARM(generate_cell)(grd_cell, r, PI,
                                                 PREC(2.0) * PI);


                        /* To get a non-symmetric grid, we simply add
                         * some more or less random number to the first
                         * latitude */
                        REAL break_symm = PREC(0.0);
                        if (s == 0)
                            break_symm = (REAL)(BREAK_SYMM);
                        grd_cell->latmax[0] -= break_symm;


                        /* Generate output file name */
                        sprintf(file, "%s/shs_%s_nx%lu_n%zu_dr%d_fft%d"
                                      "_s%d%s",
                                FOLDER, "c", nmax, i, deltar, fft,
                                (s == 0) ? 0 : 1, FTYPE);


                        f = (REAL *)malloc(grd_cell->ncell * sizeof(REAL));
                        if (f == NULL)
                        {
                            fprintf(stderr, "malloc failure.\n");
                            exit(CHARM_FAILURE);
                        }


                        CHARM(shs_cell)(grd_cell, shcs_pot, nmax, f, err);
                        CHARM(err_handler)(err, 1);
#ifdef GENREF
                        e += array2file(file, f, grd_cell->ncell);
#else
                        e += validate(file, f, grd_cell->ncell,
                                      CHARM(glob_threshold2));
#endif


                        CHARM(crd_cell_free)(grd_cell);


                        free(f);
                    }
                }
            }
        }
    }
    }
    /* ..................................................................... */






    /* Scattered cells */
    /* ..................................................................... */
    {
    size_t nlat[NSCTR] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    size_t nlon[NSCTR] = {1, 2, 3, 4, 5, 6, 7, 8, 9};


    CHARM(cell) *sctr_cell = NULL;
    for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
    {
        for (size_t i = 0; i < NSCTR; i++)
        {
            for (int deltar = 0; deltar < NDELTAR; deltar++)
            {
                REAL r = shcs_pot->r + (REAL)(DELTAR) * (REAL)deltar;


                sctr_cell = CHARM(crd_cell_calloc)(CHARM_CRD_CELL_SCATTERED,
                                                   nlat[i], nlon[i]);
                if (sctr_cell == NULL)
                {
                    fprintf(stderr, "Failed to initialize a \"crd\" "
                            "structure\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(generate_cell)(sctr_cell, r, PI, PREC(2.0) * PI);


                /* Generate output file name */
                sprintf(file, "%s/shs_%s_nx%lu_n%zu_dr%d_sctr%s",
                        FOLDER, "c", nmax, i, deltar, FTYPE);


                f = (REAL *)malloc(sctr_cell->ncell * sizeof(REAL));
                if (f == NULL)
                {
                    fprintf(stderr, "malloc failure.\n");
                    exit(CHARM_FAILURE);
                }


                CHARM(shs_cell)(sctr_cell, shcs_pot, nmax, f, err);
                CHARM(err_handler)(err, 1);
#ifdef GENREF
                e += array2file(file, f, sctr_cell->ncell);
#else
                e += validate(file, f, sctr_cell->ncell,
                              CHARM(glob_threshold2));
#endif


                CHARM(crd_cell_free)(sctr_cell);


                free(f);
            }
        }
    }
    }
    /* --------------------------------------------------------------------- */






    CHARM(err_free)(err);
    CHARM(shc_free)(shcs_pot);


    return e;
}
