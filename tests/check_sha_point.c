/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "cmp_arrays.h"
#include "parameters.h"
#include "validate.h"
/* ------------------------------------------------------------------------- */






long int check_sha_point(void)
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






    /* Check spherical harmonic analysis with points values using all supported
     * quadratures */
    /* --------------------------------------------------------------------- */
    CHARM(point) *grd_pnt = NULL;
    REAL *f               = NULL;
    CHARM(shc) *shcs_out  = NULL;
    REAL *dda             = (REAL *)malloc((SHCS_NMAX_POT + 1) * sizeof(REAL));
    if (dda == NULL)
    {
        fprintf(stderr, "Failed to initialize an array of degree "
                        "variances.\n");
        exit(CHARM_FAILURE);
    }
    REAL *dda_ref         = (REAL *)calloc((SHCS_NMAX_POT + 1), sizeof(REAL));
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


    /* Loop over all supported quadratures:
     * * "0" -- GL,
     * * "1" -- DH1,
     * * "2" -- DH2. */
    long int e = 0;
    for (int i = 0; i < 3; i++)
    {
        /* Now loop over maximum harmonic degrees */
        for (unsigned long nmax_tmp = 0; nmax_tmp <= SHCS_NMAX_POT; nmax_tmp++)
        {
            /* Create the grids on a sphere with the radius "shcs_ref->r
             * + deltar" */
            if (i == 0)
            {
                grd_pnt = CHARM(crd_point_gl)(nmax_tmp, shcs_ref->r + deltar);
                if (grd_pnt == NULL)
                {
                    fprintf(stderr, "Failed to initialize the Gauss--Legendre "
                                    "grid.\n");
                    exit(CHARM_FAILURE);
                }
            }
            else if (i == 1)
            {
                grd_pnt = CHARM(crd_point_dh1)(nmax_tmp, shcs_ref->r + deltar);
                if (grd_pnt == NULL)
                {
                    fprintf(stderr, "Failed to initialize the Driscoll--Healy "
                                    "grid (DH1).\n");
                    exit(CHARM_FAILURE);
                }
            }
            else if (i == 2)
            {
                grd_pnt = CHARM(crd_point_dh2)(nmax_tmp, shcs_ref->r + deltar);
                if (grd_pnt == NULL)
                {
                    fprintf(stderr, "Failed to initialize the Driscoll--Healy "
                                    "grid (DH2).\n");
                    exit(CHARM_FAILURE);
                }
            }


            shcs_out = CHARM(shc_calloc)(nmax_tmp, shcs_ref->mu, shcs_ref->r);
            if (shcs_out == NULL)
            {
                fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
                exit(CHARM_FAILURE);
            }


            f = (REAL *)malloc(grd_pnt->nlat * grd_pnt->nlon * sizeof(REAL));
            if (f == NULL)
            {
                fprintf(stderr, "Failed to initialize an array to store the "
                                "synthesized signal.\n");
                exit(CHARM_FAILURE);
            }


            /* We do the synthesis on the sphere with the radius "shcs_ref->r
             * + deltar" */
            CHARM(shs_point)(grd_pnt, shcs_ref, nmax_tmp, f, err);
            CHARM(err_handler)(err, 1);


            /* Now do the analysis on the sphere with the radius "shcs_out->r
             * + deltar" and rescale the coefficients to "shcs_out->r ==
             * shcs_ref->r", so that the coefficients can be validated with
             * respect to "shcs_ref" */
            CHARM(sha_point)(grd_pnt, f, nmax_tmp, shcs_out, err);
            CHARM(err_handler)(err, 1);


            CHARM(shc_dda)(shcs_ref, shcs_out, nmax_tmp, dda, err);
            CHARM(err_handler)(err, 1);


            e += cmp_arrays(dda, dda_ref, nmax_tmp + 1,
                            PREC(10.0) * CHARM(glob_threshold));


            CHARM(shc_free)(shcs_out);
            CHARM(crd_point_free)(grd_pnt);
            free(f);
        }
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    CHARM(shc_free)(shcs_ref);
    CHARM(err_free)(err);
    free(dda);
    free(dda_ref);


    return e;
    /* --------------------------------------------------------------------- */
}
