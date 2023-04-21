/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "parameters.h"
#include "validate.h"
#include "generate_cell.h"
/* ------------------------------------------------------------------------- */






long int check_sha_cell(void)
{
    /* Error structure */
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize an \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */






    /* Read the test coefficients from a text file */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_ref = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_ref == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs_ref, err);
    CHARM(err_handler)(err, 1);


    /* Modify coefficients of degrees "0" and "1" to allow for an accurate
     * validation in all precisions. */
    shcs_ref->c[0][0 - 0] = (REAL)C00;
    shcs_ref->c[0][1 - 0] = (REAL)C10;
    shcs_ref->c[1][1 - 1] = (REAL)C11;
    shcs_ref->s[1][1 - 1] = (REAL)S11;
    /* --------------------------------------------------------------------- */






    /* Check SHA with cells */
    /* --------------------------------------------------------------------- */
    CHARM(cell) *grd_cell = NULL;
    REAL *f               = NULL;
    CHARM(shc) *shcs_out  = NULL;


    /* The number of latitudes and longitudes for various grids */
    size_t nlat[NCUSTOM_GRD] = {1, 1, 3, 10};
    size_t nlon[NCUSTOM_GRD] = {1, 2, 8, 22};


    char file_c[NSTR_LONG], file_s[NSTR_LONG];


    long int e = 0;


    for (unsigned long nmax = 0; nmax <= NMAX; nmax++)
    {
        for (size_t i = 0; i < NCUSTOM_GRD; i++)
        {
            for (int s = 0; s < 2; s++) /* "s == 0" for non-symm grds */
            {
                if ((nlat[i] == 1) && (s == 1))
                    /* If there is only a single latitude, the grid is
                     * automatically non-symmetric */
                    continue;


                for (int dr = 0; dr < NDELTAR; dr++)
                {
                    REAL r = shcs_ref->r + (REAL)(DELTAR) * (REAL)(dr);


                    grd_cell = CHARM(crd_cell_calloc)(CHARM_CRD_CELL_GRID,
                                                      nlat[i], nlon[i]);
                    if (grd_cell == NULL)
                    {
                        fprintf(stderr, "Failed to initialize a "
                                        "\"crd\" structure\n");
                        exit(CHARM_FAILURE);
                    }


                    CHARM(generate_cell)(grd_cell, r, PI, PREC(2.0) * PI);


                    /* A constant to artificially get a non-symmetric grid from
                     * a symmetric grid */
                    REAL break_symm = PREC(0.0);
                    if (s == 0)
                        break_symm = (REAL)(BREAK_SYMM);


                    f = (REAL *)malloc(grd_cell->nlat * grd_cell->nlon *
                                       sizeof(REAL));
                    if (f == NULL)
                    {
                        fprintf(stderr, "malloc failure.\n");
                        exit(CHARM_FAILURE);
                    }


                    if (nlat[i] >= (nmax + 1))
                    {
                        shcs_out = CHARM(shc_calloc)(nmax, shcs_ref->mu,
                                                     shcs_ref->r);
                        if (shcs_out == NULL)
                        {
                            fprintf(stderr, "Failed to initialize a "
                                            "\"shc\" structure\n");
                            exit(CHARM_FAILURE);
                        }


                        if (s == 0)
                        {
                            if (nmax == 0)
                                continue;


                            grd_cell->latmin[0] += break_symm;
                            grd_cell->latmax[1] += break_symm;
                        }


                        CHARM(shs_cell)(grd_cell, shcs_ref, nmax, f, err);
                        CHARM(err_handler)(err, 1);


                        CHARM(sha_cell)(grd_cell, f, nmax, CHARM_SHA_CELL_AQ,
                                        shcs_out, err);
                        CHARM(err_handler)(err, 1);


                        /* Generate output file names */
                        sprintf(file_c,
                                "%s/sha_c_nx%lu_n%zu_dr%d_s%d_c%s",
                                FOLDER, nmax, i, dr,
                                (s == 0) ? 0 : 1, FTYPE);
                        sprintf(file_s,
                                "%s/sha_c_nx%lu_n%zu_dr%d_s%d_s%s",
                                FOLDER, nmax, i, dr,
                                (s == 0) ? 0 : 1, FTYPE);


                        e += validate(file_c, shcs_out->c[0],
                                      ((nmax + 2) * (nmax + 1)) / 2,
                                      PREC(10.0) * CHARM(glob_threshold2));
                        e += validate(file_s, shcs_out->s[0],
                                      ((nmax + 2) * (nmax + 1)) / 2,
                                      PREC(10.0) * CHARM(glob_threshold2));


                        CHARM(shc_free)(shcs_out);
                    }


                    CHARM(crd_cell_free)(grd_cell);
                    free(f);
                }
            }
        }
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    CHARM(shc_free)(shcs_ref);
    CHARM(err_free)(err);


    return e;
    /* --------------------------------------------------------------------- */
}
