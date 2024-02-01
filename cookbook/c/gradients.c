#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <charm/charm.h>


int main(void)
{
    /* INPUTS */
    /* ===================================================================== */
    /* Path to the input file with spherical harmonic coefficients in the "gfc"
     * format. */
    char shcs_file[] = "../../data/input/EGM96-degree10.gfc";
    /* ===================================================================== */


    /* ===================================================================== */
    /* Initialize a "charm_err" structure */
    charm_err *err = charm_err_init();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* First, let's get the maximum harmonic degree stored in "shcs_file".
     * This is especially useful to read all coefficients in "shcs_file"
     * without knowing its maximum harmonic degree a priori. */
    unsigned long nmax = charm_shc_read_gfc(shcs_file, CHARM_SHC_NMAX_MODEL,
                                            NULL, NULL, err);
    charm_err_handler(err, 1);


    /* Allocate memory for the spherical harmonic coefficients up to degree
     * "nmax" */
    charm_shc *shcs = charm_shc_malloc(nmax, 1.0, 1.0);
    if (shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Now read all coefficients in "shcs_file" */
    charm_shc_read_gfc(shcs_file, nmax, NULL, shcs, err);
    charm_err_handler(err, 1);


    /* Compute the Gauss--Legendre point grid for an "nmax" given by the
     * maximum degree stored in "shcs" and for a radius equal to the reference
     * sphere, to which the spherical harmonic coefficients are scaled to.  We
     * intentionally use here the Gauss--Legendre grid instead of the
     * Driscoll--Healy grids in order to avoid the inaccuracies due to the
     * singularities of the "y"-elements at the poles (see the
     * documentation). */
    charm_point *grd = charm_crd_point_gl(shcs->nmax, shcs->r);
    if (grd == NULL)
    {
        fprintf(stderr, "Failed to compute the Gauss--Legendre grid.\n");
        exit(CHARM_FAILURE);
    }


    /* Total number of gravitational vector elements */
    size_t np = 3;


    /* Allocate memory for the gravitational vector elements */
    double **f = (double **)malloc(np * sizeof(double *));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure\n");
        exit(CHARM_FAILURE);
    }
    for (size_t i = 0; i < np; i++)
    {
        f[i] = (double *)malloc(grd->npoint * sizeof(double));
        if (f[i] == NULL)
        {
            fprintf(stderr, "malloc failure\n");
            exit(CHARM_FAILURE);
        }
    }


    /* Compute the full first-order gradient (vector) */
    charm_shs_point_grad1(grd, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Allocate memory for the magnitude of the gravitational vector */
    double *g = (double *)malloc(grd->npoint * sizeof(double));
    if (g == NULL)
    {
        fprintf(stderr, "malloc failure\n");
        exit(CHARM_FAILURE);
    }


    /* Compute the magnitude of the gravitational acceleration (no contribution
     * due to the centrifugal force is considered here) */
    for (size_t i = 0; i < grd->npoint; i++)
        g[i] = sqrt(f[0][i] * f[0][i] + f[1][i] * f[1][i] + f[2][i] * f[2][i]);


    /* Free the memory associated with the gravitational vector */
    for (size_t i = 0; i < np; i++)
        free(f[i]);
    free(f);


    /* Allocate memory for the gravitational tensor elements */
    np = 6;  /* There are six tensor elements */
    f = (double **)malloc(np * sizeof(double *));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure\n");
        exit(CHARM_FAILURE);
    }
    for (size_t i = 0; i < np; i++)
    {
        f[i] = (double *)malloc(grd->npoint * sizeof(double));
        if (f[i] == NULL)
        {
            fprintf(stderr, "malloc failure\n");
            exit(CHARM_FAILURE);
        }
    }


    /* Now compute the full second-order gradient (tensor) */
    charm_shs_point_grad2(grd, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Check whether "fxx + fyy + fzz" is zero within numerical errors */
    double trace_error = fabs(f[0][0] + f[3][0] + f[5][0]);  /* Trace at the
                                                              * first grid
                                                              * point of "grd"
                                                              * */
    double tmp;
    for (size_t i = 1; i < grd->npoint; i++)
        tmp = fabs(f[0][i] + f[3][i] + f[5][i]);
        if (tmp > trace_error)
            trace_error = tmp;


    printf("The largest error of the gravitational tensor trace is "
           "%0.16e s^-2.\n", trace_error);


    /* Free the heap memory */
    charm_shc_free(shcs);
    charm_crd_point_free(grd);
    charm_err_free(err);
    for (size_t i = 0; i < np; i++)
        free(f[i]);
    free(f);
    free(g);


    printf("Great, all done!\n");
    /* ===================================================================== */


    exit(CHARM_SUCCESS);
}
