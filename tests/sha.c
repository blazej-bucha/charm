/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../src/prec.h"
#include "cmp_arrays.h"
#include "parameters.h"
#include "validate.h"
#include "generate_crd.h"
/* ------------------------------------------------------------------------- */






/* Tests the "sha" module */
int sha(unsigned long nmax, char SHCs_file[])
{
    /* --------------------------------------------------------------------- */
    if (nmax != 10)
    {
        fprintf(stderr, "\"nmax\" has to be \"10\".\n");
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */






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
    CHARM(shc) *shcs_ref = CHARM(shc_init)(nmax, PREC(1.0), PREC(1.0));
    if (shcs_ref == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    FILE *fptr = fopen(SHCs_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n", SHCs_file);
        exit(CHARM_FAILURE);
    }
    CHARM(shc_read_mtx)(fptr, nmax, shcs_ref, err);
    fclose(fptr);
    CHARM(err_handler)(err, 1);


    /* Modify coefficients of degrees "0" and "1" to allow for an accurate
     * validation in all precisions. */
    shcs_ref->c[0][0 - 0] = (REAL)C00;
    shcs_ref->c[0][1 - 0] = (REAL)C10;
    shcs_ref->c[1][1 - 1] = (REAL)C11;
    shcs_ref->s[1][1 - 1] = (REAL)S11;
    /* --------------------------------------------------------------------- */






    /* Check spherical harmonic analysis with points values using all supported
     * quadratures */
    /* --------------------------------------------------------------------- */
    CHARM(crd) *grd;
    REAL *f;
    CHARM(shc) *shcs_out;
    REAL *dda = (REAL *)malloc((nmax + 1) * sizeof(REAL));
    if (dda == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "variances.\n");
        exit(CHARM_FAILURE);
    }
    REAL *dda_ref = (REAL *)calloc((nmax + 1), sizeof(REAL));
    if (dda_ref == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "variances.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesize in "deltar" height above the reference sphere, then do
     * the analysis on the same sphere and rescale the coefficients to
     * "shcs_ref->r", so that the results can be compared with respect to the
     * input coefficients. */
    REAL deltar = (REAL)(DELTAR);


    /* Loop over all supported quadratures */
    int errnum = 0;
    for (int i = 0; i < 3; i++)
    {
        if (i == 0)
            printf("    Surface SHA with a Gauss--Legendre point grid...\n");
        else if (i == 1)
            printf("    Surface SHA with a Driscoll--Healy point grid "
                   "(DH1)...\n");
        else if (i == 2)
            printf("    Surface SHA with a Driscoll--Healy point grid "
                   "(DH2)...\n");


        /* Now loop over maximum harmonic degrees */
        for (unsigned long nmax_tmp = 0; nmax_tmp <= nmax; nmax_tmp++)
        {
            /* Create the grids on a sphere with the radius "shcs_ref->r
             * + deltar" */
            if (i == 0)
            {
                grd = CHARM(crd_gl)(nmax_tmp, shcs_ref->r + deltar);
                if (grd == NULL)
                {
                    fprintf(stderr, "Failed to initialize the Gauss--Legendre "
                                    "grid.\n");
                    exit(CHARM_FAILURE);
                }
            }
            else if (i == 1)
            {
                grd = CHARM(crd_dh1)(nmax_tmp, shcs_ref->r + deltar);
                if (grd == NULL)
                {
                    fprintf(stderr, "Failed to initialize the Driscoll--Healy "
                                    "grid (DH1).\n");
                    exit(CHARM_FAILURE);
                }
            }
            else if (i == 2)
            {
                grd = CHARM(crd_dh2)(nmax_tmp, shcs_ref->r + deltar);
                if (grd == NULL)
                {
                    fprintf(stderr, "Failed to initialize the Driscoll--Healy "
                                    "grid (DH2).\n");
                    exit(CHARM_FAILURE);
                }
            }


            shcs_out = CHARM(shc_init)(nmax_tmp, shcs_ref->mu, shcs_ref->r);
            if (shcs_out == NULL)
            {
                fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
                exit(CHARM_FAILURE);
            }


            f = (REAL *)malloc(grd->nlat * grd->nlon * sizeof(REAL));
            if (f == NULL)
            {
                fprintf(stderr, "Failed to initialize an array to store the "
                                "synthesized signal.\n");
                exit(CHARM_FAILURE);
            }


            /* We do the synthesis on the sphere with the radius "shcs_ref->r
             * + deltar" */
            CHARM(shs_point)(grd, shcs_ref, nmax_tmp, f, err);
            CHARM(err_handler)(err, 1);


            /* Now do the analysis on the sphere with the radius "shcs_out->r
             * + deltar" and rescale the coefficients to "shcs_out->r ==
             * shcs_ref->r", so that the coefficients can be validated with
             * respect to "shcs_ref" */
            CHARM(sha_point)(grd, f, nmax_tmp, shcs_out, err);
            CHARM(err_handler)(err, 1);


            CHARM(shc_dda)(shcs_ref, shcs_out, nmax_tmp, dda, err);
            CHARM(err_handler)(err, 1);


            errnum += cmp_arrays(dda, dda_ref, nmax_tmp + 1,
                                 PREC(10.0) * CHARM(glob_threshold));


            CHARM(shc_free)(shcs_out);
            CHARM(crd_free)(grd);
            free(f);
        }
    }


    free(dda);
    free(dda_ref);
    /* --------------------------------------------------------------------- */






    /* Check SHA with cells */
    /* --------------------------------------------------------------------- */
    /* The number of latitudes and longitudes for various grids */
    size_t nlat[NCUSTOM_GRD] = {1, 1, 3, 10};
    size_t nlon[NCUSTOM_GRD] = {1, 2, 8, 22};


    char file_c[NSTR], file_s[NSTR];


    printf("    Surface SHA with cell grids...\n");
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


                    CHARM(crd) *grd = CHARM(crd_init)(CHARM_CRD_CELLS_GRID,
                                                      nlat[i], nlon[i]);
                    if (grd == NULL)
                    {
                        fprintf(stderr, "Failed to initialize a "
                                        "\"crd\" structure\n");
                        exit(1);
                    }


                    CHARM(generate_crd)(grd, r, PI, PREC(2.0) * PI);


                    /* A constant to artificially get a non-symmetric grid from
                     * a symmetric grid */
                    REAL break_symm = PREC(0.0);
                    if (s == 0)
                        break_symm = (REAL)(BREAK_SYMM);


                    REAL *f = (REAL *)malloc(grd->nlat * grd->nlon *
                                             sizeof(REAL));
                    if (f == NULL)
                    {
                        fprintf(stderr, "malloc failure.\n");
                        exit(1);
                    }


                    if (nlat[i] >= (nmax + 1))
                    {
                        shcs_out = CHARM(shc_init)(nmax, shcs_ref->mu,
                                                         shcs_ref->r);
                        if (shcs_out == NULL)
                        {
                            fprintf(stderr, "Failed to initialize a "
                                            "\"shc\" structure\n");
                            exit(1);
                        }


                        if (s == 0)
                        {
                            if (nmax == 0)
                                continue;


                            grd->lat[1] += break_symm;
                            grd->lat[2] += break_symm;
                        }


                        CHARM(shs_cell)(grd, shcs_ref, nmax, f, err);
                        CHARM(err_handler)(err, 1);


                        CHARM(sha_cell)(grd, f, nmax, CHARM_SHA_CELL_AQ,
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


                        errnum += validate(file_c, shcs_out->c[0],
                                           ((nmax + 2) * (nmax + 1)) / 2,
                                           PREC(10.0) *
                                           CHARM(glob_threshold2));
                        errnum += validate(file_s, shcs_out->s[0],
                                           ((nmax + 2) * (nmax + 1)) / 2,
                                           PREC(10.0) *
                                           CHARM(glob_threshold2));


                        CHARM(shc_free)(shcs_out);
                    }


                    CHARM(crd_free)(grd);
                    free(f);
                }
            }
        }
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    CHARM(shc_free)(shcs_ref);
    CHARM(err_free)(err);


    return errnum;
    /* --------------------------------------------------------------------- */
}

