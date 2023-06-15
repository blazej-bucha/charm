#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <charm/charm.h>


int main(void)
{
    /* INPUTS */
    /* ===================================================================== */
    /* Define the path to an input text file with spherical harmonic
     * coefficients.  For details on the structure of the text file, see the
     * description of the "charm_shc_read_mtx" function in the "charm_shc"
     * module. */
    char shcs_in_file[] = "../../data/input/EGM96-degree10-mtx.txt";


    /* Maximum harmonic degree of coefficients to read.  The same degree is
     * used later in the harmonic synthesis and harmonic analysis. */
    unsigned long nmax = 10;
    /* ===================================================================== */






    /* ===================================================================== */
    printf("Closed-loop experiment -- Spherical harmonic synthesis and "
           "analysis\n");
    printf("===========================\n");


    /* Initialize a "charm_shc" structure */
    charm_shc *shcs = charm_shc_calloc(nmax, 1.0, 1.0);
    if (shcs == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Initialize a "charm_err" structure */
    charm_err *err = charm_err_init();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Read spherical harmonic coefficients from the input text file. */
    charm_shc_read_mtx(shcs_in_file, nmax, shcs, err);
    charm_err_handler(err, 1);


    /* Now let's say we do not want to use the zero-degree term.  This can be
     * achieved simply by setting the respective "C00" coefficient to zero. */
    shcs->c[0][0] = 0.0;


    /* Compute the Gauss--Legendre grid for a given "nmax" on a sphere that is
     * "1000.0" metres above the reference sphere of spherical harmonic
     * coefficients. */
    charm_point *grd_pnt = charm_crd_point_gl(nmax, shcs->r + 1000.0);
    if (grd_pnt == NULL)
    {
        fprintf(stderr, "Failed to compute the Gauss--Legendre grid.\n");
        exit(CHARM_FAILURE);
    }


    /* Allocate the memory for the synthesized signal */
    double *f = (double *)malloc(grd_pnt->nlat * grd_pnt->nlon *
                                 sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Perform the synthesis.  Since the Gauss--Legendre grid resides "1000.0"
     * metres above the reference sphere of the coefficients, performed is
     * solid spherical harmonic synthesis. */
    charm_shs_point(grd_pnt, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Print some synthesized values */
    printf("Print some synthesized values of the signal...\n");
    size_t i = 0;
    size_t j = 0;
    printf("f(%zu, %zu) = %0.16e\n", i, j, f[i * grd_pnt->nlon + j]);


    i = grd_pnt->nlat / 2;
    j = grd_pnt->nlon / 2;
    printf("f(%zu, %zu) = %0.16e\n", i, j, f[i * grd_pnt->nlon + j]);


    /* Initialize a new structure of spherical harmonic coefficients for
     * coefficients to be computed by the harmonic analysis.  Used are the same
     * scalling constants as in "shcs". */
    charm_shc *shcs2 = charm_shc_calloc(nmax, shcs->mu, shcs->r);
    if (shcs2 == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Now use the synthesized signal and compute back its spherical harmonic
     * coefficients by harmonic analysis.  The output coefficients in "shcs2"
     * should be the same as the input ones "shcs".  Note that the signal "f"
     * resides on a sphere that is "1000.0" metres above the reference sphere.
     * Therefore, the coefficients are automatically properly rescaled to the
     * desired sphere defined by "shcs2->r". */
    charm_sha_point(grd_pnt, f, nmax, shcs2, err);
    charm_err_handler(err, 1);


    /* Now check whether "shcs" and "shcs2" are the same by computing their
     * difference degree amplitudes */
    double *dda = (double *)malloc((nmax + 1) * sizeof(double));
    if (dda == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }
    charm_shc_dda(shcs, shcs2, nmax, dda, err);
    charm_err_handler(err, 1);


    /* Print some difference degree amplitudes */
    unsigned long n = 0;
    printf("\n\n");
    printf("Now print the difference degree amplitudes. ");
    printf("These should be very small, say, at the order of 1e-18 "
           "or less...\n");
    printf("Difference degree amplitude for harmonic degree %lu = %0.16e\n",
            n, dda[n]);
    n = 4;
    printf("Difference degree amplitude for harmonic degree %lu = %0.16e\n",
            n, dda[n]);
    n = 10;
    printf("Difference degree amplitude for harmonic degree %lu = %0.16e\n",
            n, dda[n]);


    /* We will not need some of the structures and arrays anymore, so let's
     * free them. */
    charm_crd_point_free(grd_pnt);
    charm_shc_free(shcs2);
    free(f);
    free(dda);


    printf("===========================\n\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Solid spherical harmonic synthesis at scattered points\n");
    printf("===========================\n");


    /* Let's create a "charm_point" structure to store 3 scattered points. */
    charm_point *sctr_pnt = charm_crd_point_calloc(CHARM_CRD_POINT_SCATTERED,
                                                   3, 3);
    if (sctr_pnt == NULL)
    {
        fprintf(stderr,
                "Failed to initialize the \"charm_point\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Now let's feed "sctr_pnt" with some computation points (angular values
     * must be provided in radians). */
    sctr_pnt->lat[0] =  0.1;
    sctr_pnt->lat[1] =  0.436231;
    sctr_pnt->lat[2] = -0.9651;
    sctr_pnt->lon[0] =  0.0;
    sctr_pnt->lon[1] =  1.53434;
    sctr_pnt->lon[2] =  4.2316;
    sctr_pnt->r[0]   = shcs->r + 1000.0;
    sctr_pnt->r[1]   = shcs->r + 2000.0;
    sctr_pnt->r[2]   = shcs->r + 3000.0;


    /* Allocate memory to store the synthesized signal */
    f = (double *)malloc(sctr_pnt->nlat * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis */
    charm_shs_point(sctr_pnt, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* It's really this easy! */


    /* Print the synthesized values */
    for (size_t i = 0; i < sctr_pnt->nlat; i++)
        printf("%0.16e\n", f[i]);


    charm_crd_point_free(sctr_pnt);
    free(f);


    printf("===========================\n\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Solid spherical harmonic synthesis at a custom grid of points\n");
    printf("===========================\n");


    /* Initialize a charm_point structure to hold a "5 x 10" point grid. */
    grd_pnt = charm_crd_point_calloc(CHARM_CRD_POINT_GRID, 5, 10);
    if (grd_pnt == NULL)
    {
        fprintf(stderr,
                "Failed to initialize the \"charm_point\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Define some grid points */
    for (size_t i = 0; i < grd_pnt->nlat; i++)
    {
        grd_pnt->lat[i] = M_PI_2 - ((double)i / (double)grd_pnt->nlat) * M_PI;
        grd_pnt->r[i]   = shcs->r + (double)i;
    }
    for (size_t j = 0; j < grd_pnt->nlon; j++)
    {
        grd_pnt->lon[j] = ((double)j / (double)grd_pnt->nlon) * (2.0 * M_PI);
    }


    /* Initialize an array to store the synthesized signal */
    f = (double *)malloc(grd_pnt->nlat * grd_pnt->nlon * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis */
    charm_shs_point(grd_pnt, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Print the synthesized values */
    for (size_t i = 0; i < grd_pnt->nlat; i++)
        for (size_t j = 0; j < grd_pnt->nlon; j++)
            printf("%0.16e\n", f[i * grd_pnt->nlon + j]);


    charm_crd_point_free(grd_pnt);
    free(f);


    printf("===========================\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Solid spherical harmonic synthesis at scattered cells\n");
    printf("===========================\n");


    /* Initialize a "charm_cell" structure to hold 3 scattered cells. */
    charm_cell *sctr_cell = charm_crd_cell_calloc(CHARM_CRD_CELL_SCATTERED,
                                                  3, 3);
    if (sctr_cell == NULL)
    {
        fprintf(stderr, "Failed to initialize scattered cells.\n");
        exit(CHARM_FAILURE);
    }


    /* Define more or less randomly some scattered cells */
    /* --------------------------------------------------------------------- */
    /* The maximum latitude of the first cell */
    sctr_cell->latmax[0] =  0.323413435;

    /* The minimum latitude of the first cell */
    sctr_cell->latmin[0] =  sctr_cell->latmax[0] - 0.234323;

    /* The maximum latitude of the second cell */
    sctr_cell->latmax[1] = -0.90234320952;

    /* The minimum latitude of the second cell */
    sctr_cell->latmin[1] =  sctr_cell->latmax[1] - 0.4456;

    /* And so on... */
    sctr_cell->latmax[2] =  0.0;
    sctr_cell->latmin[2] =  sctr_cell->latmax[2] - M_PI_2;

    /* And now the same, but with the longitudes... */
    sctr_cell->lonmin[0] =  0.123456789;
    sctr_cell->lonmax[0] =  sctr_cell->lonmin[0] + 1.3235;
    sctr_cell->lonmin[1] =  4.3445234;
    sctr_cell->lonmax[1] =  sctr_cell->lonmin[1] + 0.01;
    sctr_cell->lonmin[2] =  0.0;
    sctr_cell->lonmax[2] =  sctr_cell->lonmin[2] + 2.0 * M_PI;

    /* Finally, we define spherical radii of the scattered cells.  Importantly,
     * the spherical radius is *constants* over the cells (but may vary from
     * cell to cell). */
    sctr_cell->r[0]      =  shcs->r;
    sctr_cell->r[1]      =  shcs->r + 1000.0;
    sctr_cell->r[2]      =  shcs->r + 2000.0;
    /* --------------------------------------------------------------------- */


    /* Initialize an array to store the synthesized signal */
    f = (double *)malloc(sctr_cell->nlat * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis */
    charm_shs_cell(sctr_cell, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Print the synthesized values */
    for (size_t i = 0; i < sctr_cell->nlat; i++)
        printf("%0.16e\n", f[i]);


    charm_crd_cell_free(sctr_cell);
    free(f);


    printf("===========================\n\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Solid spherical harmonic synthesis at a grid "
           "of cells\n");
    printf("===========================\n");


    /* Initialize a "charm_cell" structure to hold a grid of cells with "15"
     * cells in the latitudinal direction and "30" cells in the longitudinal
     * direction. */
    charm_cell *grd_cell = charm_crd_cell_calloc(CHARM_CRD_CELL_GRID, 15, 30);


    /* Define some grid cells */
    for (size_t i = 0; i < grd_cell->nlat; i++)
    {
        /* The maximum cell latitudes */
        grd_cell->latmax[i] = M_PI_2 -
                              ((double)i / (double)grd_cell->nlat) * M_PI;

        /* The minimum cell latitudes */
        grd_cell->latmin[i] = M_PI_2 -
                              ((double)(i + 1) / (double)grd_cell->nlat) *
                              M_PI;

        /* The spherical radii (may vary with the latitude, but we use
         * a constant value here, because of the example that follows after
         * this one). */
        grd_cell->r[i] = shcs->r + 1000.0;
    }
    for (size_t j = 0; j < grd_cell->nlon; j++)
    {
        /* The minimum cell longitudes */
        grd_cell->lonmin[j] = ((double)j / (double)grd_cell->nlon) *
                              (2.0 * M_PI);

        /* The maximum cell longitudes */
        grd_cell->lonmax[j] = ((double)(j + 1) / (double)grd_cell->nlon) *
                              (2.0 * M_PI);
    }


    /* Initialize an array to store the synthesized signal */
    f = (double *)malloc(grd_cell->nlat * grd_cell->nlon * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis */
    charm_shs_cell(grd_cell, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Print some of the synthesized values */
    /* --------------------------------------------------------------------- */
    i = 0;
    j = 10;
    printf("f(%zu, %zu) = %0.16e\n", i, j, f[i * grd_cell->nlon + j]);


    i = 4;
    j = 3;
    printf("f(%zu, %zu) = %0.16e\n", i, j, f[i * grd_cell->nlon + j]);
    /* --------------------------------------------------------------------- */


    printf("===========================\n\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Spherical harmonic analysis of block-mean values in cells\n");
    printf("===========================\n");


    /* In this example, we use the "grd_cell" structure and the signal "f" from
     * the previous example. */


    /* Initialize a new structure of spherical harmonic coefficients for
     * coefficients to be recovered from the harmonic analysis.  Used are the
     * same scalling constants as in "shcs". */
    shcs2 = charm_shc_calloc(nmax, shcs->mu, shcs->r);
    if (shcs2 == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the analysis using the method of the approximate quadrature */
    charm_sha_cell(grd_cell, f, nmax, CHARM_SHA_CELL_AQ, shcs2, err);
    charm_err_handler(err, 1);


    /* Print some of the computed coefficients.  Note that the harmonic
     * analysis with block-mean values in cells in *not* exact, hence the
     * coefficients will not be equal to the input ones. */
    /* --------------------------------------------------------------------- */
    n = 2;
    unsigned long m = 0;
    printf("C(%3lu,%3lu) = %0.16e\n", n, m, shcs2->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.16e\n", n, m, shcs2->s[m][n - m]);


    n = 9;
    m = 0;
    printf("C(%3lu,%3lu) = %0.16e\n", n, m, shcs2->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.16e\n", n, m, shcs2->s[m][n - m]);


    n = 9;
    m = 4;
    printf("C(%3lu,%3lu) = %0.16e\n", n, m, shcs2->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.16e\n", n, m, shcs2->s[m][n - m]);


    n = 9;
    m = 9;
    printf("C(%3lu,%3lu) = %0.16e\n", n, m, shcs2->c[m][n - m]);
    printf("S(%3lu,%3lu) = %0.16e\n", n, m, shcs2->s[m][n - m]);
    /* --------------------------------------------------------------------- */


    charm_crd_cell_free(grd_cell);
    charm_shc_free(shcs2);
    free(f);


    printf("===========================\n\n");
    /* ===================================================================== */






    /* Free the heap memory */
    /* ===================================================================== */
    charm_err_free(err);
    charm_shc_free(shcs);


    printf("\nGreat, all done!\n");


    exit(CHARM_SUCCESS);
    /* ===================================================================== */
}
