/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "parameters.h"
#include "cmp_arrays.h"
/* ------------------------------------------------------------------------- */






/* Checks the routine to rescale spherical harmonic coefficients by
 *
 * 1) loading a reference set of coefficients,
 *
 * 2) synthesizing the signal from the reference coefficients,
 *
 * 3) rescaling the coefficients,
 *
 * 4) synthesizing the signal from the rescaled coefficients, and
 *
 * 5) comparing the two signals.  */
long int check_shc_rescale(void)
{
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        printf("Failed to initialize the \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shc) *shcs = CHARM(shc_calloc)(SHCS_NMAX_POT, PREC(1.0), PREC(1.0));
    if (shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize a \"shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shc_read_mtx)(SHCS_IN_PATH_POT_MTX, SHCS_NMAX_POT, shcs, err);
    CHARM(err_handler)(err, 1);


    CHARM(point) *grd = CHARM(crd_point_gl)(SHCS_NMAX_POT,
                                            shcs->r + (REAL)DELTAR);
    if (grd == NULL)
    {
        fprintf(stderr, "Failed to compute the Gauss--Legendre grid.\n");
        exit(CHARM_FAILURE);
    }


    REAL *f1 = (REAL *)malloc((grd->nlat * grd->nlon) * sizeof(REAL));
    if (f1 == NULL)
    {
        fprintf(stderr, "Failed to allocate the input signal.\n");
        exit(CHARM_FAILURE);
    }


    REAL *f2 = (REAL *)malloc((grd->nlat * grd->nlon) * sizeof(REAL));
    if (f2 == NULL)
    {
        fprintf(stderr, "Failed to allocate the input signal.\n");
        exit(CHARM_FAILURE);
    }


    CHARM(shs_point)(grd, shcs, SHCS_NMAX_POT, f1, err);
    CHARM(err_handler)(err, 1);


    CHARM(shc_rescale)(shcs, shcs->mu * PREC(1.1), shcs->r * PREC(1.1), err);
    CHARM(err_handler)(err, 1);


    CHARM(shs_point)(grd, shcs, SHCS_NMAX_POT, f2, err);
    CHARM(err_handler)(err, 1);


    long int e = cmp_arrays(f1, f2, grd->nlat * grd->nlon,
                            PREC(10.0) * CHARM(glob_threshold));


    CHARM(crd_point_free)(grd);
    CHARM(err_free)(err);
    CHARM(shc_free)(shcs);
    free(f1);
    free(f2);


    return e;
}
