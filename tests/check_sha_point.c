/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "cmp_arrays.h"
#include "parameters.h"
#include "error_messages.h"
#include "modify_low_degree_coefficients.h"
#include "check_sha_point.h"
/* ------------------------------------------------------------------------- */






long int check_sha_point(void)
{
    /* Error structure */
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */






    /* Read the test coefficients from a text file */
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs_ref = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0),
                                             PREC(1.0));
    if (shcs_ref == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs_ref, err);
    CHARM(err_handler)(err, 1);


    /* Modify coefficients of degrees "0" and "1" to allow for an accurate
     * validation in all precisions. */
    modify_low_degree_coefficients(shcs_ref);
    /* --------------------------------------------------------------------- */






    /* Check spherical harmonic analysis with points values using all supported
     * quadratures */
    /* --------------------------------------------------------------------- */
    CHARM(point) *grd_pnt = NULL;
    REAL *f               = NULL;
    CHARM(shc) *shcs_out  = NULL;
    REAL *dda = (REAL *)malloc((SHCS_NMAX_POT + 1) * sizeof(REAL));
    if (dda == NULL)
    {
        fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }
    REAL *dda_ref = (REAL *)calloc((SHCS_NMAX_POT + 1), sizeof(REAL));
    if (dda_ref == NULL)
    {
        fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
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
                grd_pnt = CHARM(crd_point_gl)(nmax_tmp, shcs_ref->r + deltar);
            else if (i == 1)
                grd_pnt = CHARM(crd_point_dh1)(nmax_tmp, shcs_ref->r + deltar);
            else if (i == 2)
                grd_pnt = CHARM(crd_point_dh2)(nmax_tmp, shcs_ref->r + deltar);
            if (grd_pnt == NULL)
            {
                fprintf(stderr, "%s", ERR_MSG_POINT);
                exit(CHARM_FAILURE);
            }


            shcs_out = CHARM(shc_calloc)(nmax_tmp, shcs_ref->mu, shcs_ref->r);
            if (shcs_out == NULL)
            {
                fprintf(stderr, "%s", ERR_MSG_SHC);
                exit(CHARM_FAILURE);
            }


            f = (REAL *)malloc(grd_pnt->npoint * sizeof(REAL));
            if (f == NULL)
            {
                fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
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


    CHARM(shc_free)(shcs_ref);
    free(dda_ref);
    free(dda);
    /* --------------------------------------------------------------------- */






    /* Now check the dynamical switching and the loop unrolling.  To do such
     * a test, we need sufficiently large harmonic degrees, so that CHarm will
     * actually apply the dynamical switching and loop unrolling.  But to avoid
     * having large input data files, we simply set spherical harmonic
     * coefficients to some fake values in this test.  Except for this, the
     * test is almost identical to the previous ones that used actual spherical
     * harmonic coefficients. */
    /* --------------------------------------------------------------------- */
    for (int i = 0; i < 3; i++)
    {
        for (unsigned long nmax_tmp = NMAX_DS_MIN; nmax_tmp <= NMAX_DS_MAX;
             nmax_tmp++)
        {
            shcs_ref = CHARM(shc_calloc)(nmax_tmp, PREC(1.0), PREC(1.0));
            if (shcs_ref == NULL)
            {
                fprintf(stderr, "%s", ERR_MSG_SHC);
                exit(CHARM_FAILURE);
            }
            for (size_t m = 0; m <= nmax_tmp; m++)
            {
                for (size_t n = m; n <= nmax_tmp; n++)
                {
                    shcs_ref->c[m][n - m] = PREC(1.0);
                    if (m > 0)
                        shcs_ref->s[m][n - m] = PREC(1.0);
                }
            }


            dda_ref = (REAL *)calloc((nmax_tmp + 1), sizeof(REAL));
            if (dda_ref == NULL)
            {
                fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
                exit(CHARM_FAILURE);
            }


            dda = (REAL *)calloc((nmax_tmp + 1), sizeof(REAL));
            if (dda_ref == NULL)
            {
                fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
                exit(CHARM_FAILURE);
            }


            if (i == 0)
                grd_pnt = CHARM(crd_point_gl)(nmax_tmp, shcs_ref->r);
            else if (i == 1)
                grd_pnt = CHARM(crd_point_dh1)(nmax_tmp, shcs_ref->r);
            else if (i == 2)
                grd_pnt = CHARM(crd_point_dh2)(nmax_tmp, shcs_ref->r);
            if (grd_pnt == NULL)
            {
                fprintf(stderr, "%s", ERR_MSG_POINT);
                exit(CHARM_FAILURE);
            }


            shcs_out = CHARM(shc_calloc)(nmax_tmp, shcs_ref->mu, shcs_ref->r);
            if (shcs_out == NULL)
            {
                fprintf(stderr, "%s", ERR_MSG_SHC);
                exit(CHARM_FAILURE);
            }


            f = (REAL *)malloc(grd_pnt->npoint * sizeof(REAL));
            if (f == NULL)
            {
                fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
                exit(CHARM_FAILURE);
            }


            CHARM(shs_point)(grd_pnt, shcs_ref, nmax_tmp, f, err);
            CHARM(err_handler)(err, 1);


            CHARM(sha_point)(grd_pnt, f, nmax_tmp, shcs_out, err);
            CHARM(err_handler)(err, 1);


            CHARM(shc_dda)(shcs_ref, shcs_out, nmax_tmp, dda, err);
            CHARM(err_handler)(err, 1);


            /* A more relaxed threshold is needed here, given the high-degree
             * analysis and synthesis */
            e += cmp_arrays(dda, dda_ref, nmax_tmp + 1,
                            CHARM(glob_threshold2));


            CHARM(shc_free)(shcs_ref);
            CHARM(shc_free)(shcs_out);
            CHARM(crd_point_free)(grd_pnt);
            free(f);
            free(dda);
            free(dda_ref);
        }
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);


    return e;
    /* --------------------------------------------------------------------- */
}
