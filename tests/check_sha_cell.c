/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "parameters.h"
#include "error_messages.h"
#include "generate_cell.h"
#ifdef GENREF
#   include "array2file.h"
#else
#   include "validate.h"
#endif
#include "modify_low_degree_coefficients.h"
#include "check_sha_cell.h"
/* ------------------------------------------------------------------------- */






long int check_sha_cell(void)
{
    /* Error structure */
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */






    /* Read the test coefficients from a text file */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_ref = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_ref == NULL)
    {
        fprintf(stderr, ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs_ref, err);
    CHARM(err_handler)(err, 1);


    /* Modify coefficients of degrees "0" and "1" to allow for an accurate
     * validation in all precisions. */
    modify_low_degree_coefficients(shcs_ref);
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
                else if ((nlat[i] >= (nmax + 1)) && (s == 0) && (nmax == 0))
                    continue;


                for (int dr = 0; dr < NDELTAR; dr++)
                {
                    if (nlat[i] >= (nmax + 1))
                    {
                        REAL r = shcs_ref->r + (REAL)(DELTAR) * (REAL)(dr);


                        grd_cell = CHARM(crd_cell_calloc)(CHARM_CRD_CELL_GRID,
                                                          nlat[i], nlon[i]);
                        if (grd_cell == NULL)
                        {
                            fprintf(stderr, ERR_MSG_CELL);
                            exit(CHARM_FAILURE);
                        }


                        CHARM(generate_cell)(grd_cell, r, PI, PREC(2.0) * PI);


                        /* A constant to artificially get a non-symmetric grid
                         * from a symmetric grid */
                        REAL break_symm = PREC(0.0);
                        if (s == 0)
                            break_symm = (REAL)(BREAK_SYMM);


                        f = (REAL *)malloc(grd_cell->ncell * sizeof(REAL));
                        if (f == NULL)
                        {
                            fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
                            exit(CHARM_FAILURE);
                        }


                        shcs_out = CHARM(shc_calloc)(nmax, shcs_ref->mu,
                                                     shcs_ref->r);
                        if (shcs_out == NULL)
                        {
                            fprintf(stderr, ERR_MSG_SHC);
                            exit(CHARM_FAILURE);
                        }


                        if (s == 0)
                        {
                            grd_cell->latmin[0] += break_symm;
                            grd_cell->latmax[1] += break_symm;
                        }


                        CHARM(shs_cell)(grd_cell, shcs_ref, nmax, f, err);
                        CHARM(err_handler)(err, 1);


                        CHARM(sha_cell)(grd_cell, f, nmax, CHARM_SHA_CELL_AQ,
                                        shcs_out, err);
                        CHARM(err_handler)(err, 1);


                        /* Generate output file names */
                        snprintf(file_c, NSTR_LONG,
                                 "%s/sha_c_nx%lu_n%zu_dr%d_s%d_c%s",
                                 FOLDER, nmax, i, dr,
                                 (s == 0) ? 0 : 1, FTYPE);
                        snprintf(file_s, NSTR_LONG,
                                 "%s/sha_c_nx%lu_n%zu_dr%d_s%d_s%s",
                                 FOLDER, nmax, i, dr,
                                 (s == 0) ? 0 : 1, FTYPE);


#ifdef GENREF
                        e += array2file(file_c, shcs_out->c[0],
                                        ((nmax + 2) * (nmax + 1)) / 2);
                        e += array2file(file_s, shcs_out->s[0],
                                        ((nmax + 2) * (nmax + 1)) / 2);
#else
                        e += validate(file_c, shcs_out->c[0],
                                      ((nmax + 2) * (nmax + 1)) / 2,
                                      PREC(10.0) * CHARM(glob_threshold2));
                        e += validate(file_s, shcs_out->s[0],
                                      ((nmax + 2) * (nmax + 1)) / 2,
                                      PREC(10.0) * CHARM(glob_threshold2));
#endif


                        CHARM(shc_free)(shcs_out);
                        CHARM(crd_cell_free)(grd_cell);
                        free(f);
                    }
                }
            }
        }
    }


    /* Check that analysis with zero cells in "charm_cell" do not produce any
     * kind of error */
    /* ..................................................................... */
    f = NULL;
    grd_cell = CHARM(crd_cell_malloc)(CHARM_CRD_CELL_GRID, 0, 0);
    shcs_out = CHARM(shc_malloc)(SHCS_NMAX_POT, PREC(1.0), PREC(1.0));


    CHARM(sha_cell)(grd_cell, f, shcs_out->nmax, CHARM_SHA_CELL_AQ, shcs_out,
                    err);
    if (!CHARM(err_isempty)(err))
    {
        printf("\n        WARNING: Analysis with zero cells didn't pass!\n");
        e += 1;
    }


    CHARM(err_reset)(err);
    CHARM(shc_free)(shcs_out);
    CHARM(crd_cell_free)(grd_cell);
    free(f);
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    CHARM(shc_free)(shcs_ref);
    CHARM(err_free)(err);


    return e;
    /* --------------------------------------------------------------------- */
}
