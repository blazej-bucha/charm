#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
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
    char shcs_in_file[] = "../data/input/EGM96-degree10-mtx.txt";


    /* Maximum harmonic coefficients to be read.  The same degree is used later
     * in the harmonic synthesis and harmonic analysis. */
    unsigned long nmax = 10;
    /* ===================================================================== */






    /* ===================================================================== */
    printf("Closed-loop experiment -- Spherical harmonic synthesis and "
           "analysis\n");
    printf("===========================\n");


    /* Initialize a "charm_shc" structure */
    charm_shc *shcs = charm_shc_init(nmax, 1.0, 1.0);
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
    FILE *fptr = fopen(shcs_in_file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open the stream for \"%s\".\n",
                shcs_in_file);
        exit(CHARM_FAILURE);
    }
    charm_shc_read_mtx(fptr, nmax, shcs, err);
    charm_err_handler(err, 1);
    fclose(fptr);


    /* Now let's say we do not want to use the zero-degree term.  This can be
     * achieved simply by setting the respective "C00" coefficient to zero. */
    shcs->c[0][0] = 0.0;


    /* Compute the Gauss--Legendre grid for a given "nmax" on a sphere that is
     * "1000.0" metres above the reference sphere of spherical harmonic
     * coefficients. */
    charm_crd *grd = charm_crd_gl(nmax, shcs->r + 1000.0);
    if (grd == NULL)
    {
        fprintf(stderr, "Failed to compute the Gauss--Legendre grid.\n");
        exit(CHARM_FAILURE);
    }


    /* Allocate the memory for the synthesized signal */
    double *f = (double *)malloc(grd->nlat * grd->nlon * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Perform the synthesis.  Since the Gauss--Legendre grid resides "1000.0"
     * metres above the reference sphere of the coefficients, performed is
     * solid spherical harmonic synthesis. */
    charm_shs_point(grd, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Print some synthesized values */
    printf("Print some synthesized values of the signal...\n");
    size_t i = 0;
    size_t j = 0;
    printf("f(%zu, %zu) = %0.16e\n", i, j, f[i * grd->nlon + j]);


    i = grd->nlat / 2;
    j = grd->nlon / 2;
    printf("f(%zu, %zu) = %0.16e\n", i, j, f[i * grd->nlon + j]);


    /* Initialize a new structure of spherical harmonic coefficients for
     * coefficients to be computed by the harmonic analysis.  Used are the same
     * scalling constants as in "shcs". */
    charm_shc *shcs2 = charm_shc_init(nmax, shcs->mu, shcs->r);
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
    charm_sha_point(grd, f, nmax, shcs2, err);
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
    charm_crd_free(grd);
    charm_shc_free(shcs2);
    free(f);
    free(dda);


    printf("===========================\n\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Solid spherical harmonic synthesis at scattered points\n");
    printf("===========================\n");


    /* Let's create a "charm_crd" structure to store 3 scattered points. */
    charm_crd *sctr = charm_crd_init(CHARM_CRD_POINTS_SCATTERED, 3, 3);
    if (sctr == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_crd\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Now let's feed "sctr" with some computation points (angular values must
     * be provided in radians). */
    sctr->lat[0] =  0.1;
    sctr->lat[1] =  0.436231;
    sctr->lat[2] = -0.9651;
    sctr->lon[0] =  0.0;
    sctr->lon[1] =  1.53434;
    sctr->lon[2] =  4.2316;
    sctr->r[0]   = shcs->r + 1000.0;
    sctr->r[1]   = shcs->r + 2000.0;
    sctr->r[2]   = shcs->r + 3000.0;


    /* Allocate memory to store the synthesized signal */
    f = (double *)malloc(sctr->nlat * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis */
    charm_shs_point(sctr, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* It's really this easy! */


    /* Print the synthesized values */
    for (size_t i = 0; i < sctr->nlat; i++)
        printf("%0.16e\n", f[i]);


    charm_crd_free(sctr);
    free(f);


    printf("===========================\n\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Solid spherical harmonic synthesis at a custom grid of points\n");
    printf("===========================\n");


    /* Initialize a charm_crd structure to hold a "5 x 10" points grid. */
    grd = charm_crd_init(CHARM_CRD_POINTS_GRID, 5, 10);
    if (grd == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_crd\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Define some grid points */
    for (size_t i = 0; i < grd->nlat; i++)
    {
        grd->lat[i] = -M_PI_2 + (double)i * M_PI / grd->nlat;
        grd->r[i]   = shcs->r + (double)i;
    }
    for (size_t j = 0; j < grd->nlon; j++)
    {
        grd->lon[j] = (double)j * (2.0 * M_PI) / grd->nlon;
    }


    /* Initialize an array to store the synthesized signal */
    f = (double *)malloc(grd->nlat * grd->nlon * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis */
    charm_shs_point(grd, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Print the synthesized values */
    for (size_t i = 0; i < grd->nlat; i++)
        for (size_t j = 0; j < grd->nlon; j++)
            printf("%0.16e\n", f[i * grd->nlon + j]);


    charm_crd_free(grd);
    free(f);


    printf("===========================\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Solid spherical harmonic synthesis at scattered cells\n");
    printf("===========================\n");


    /* Initialize a "charm_crd" structure to hold 3 scattered cells. */
    sctr = charm_crd_init(CHARM_CRD_CELLS_SCATTERED, 3, 3);
    if (sctr == NULL)
    {
        fprintf(stderr, "Failed to initialize scattered cells.\n");
        exit(CHARM_FAILURE);
    }


    /* Define more or less randomly some scattered cells */
    /* --------------------------------------------------------------------- */
    /* The minimum latitude of the first cell */
    sctr->lat[0] =  0.323413435;

    /* The maximum latitude of the first cell */
    sctr->lat[1] =  sctr->lat[0] + 0.234323;

    /* The minimum latitude of the second cell */
    sctr->lat[2] = -0.90234320952;

    /* The maximum latitude of the second cell */
    sctr->lat[3] =  sctr->lat[2] + 0.4456;

    /* And so on... */
    sctr->lat[4] =  0.0;
    sctr->lat[5] =  sctr->lat[4] + M_PI_2;

    /* And now the same, but with the longitudes... */
    sctr->lon[0] =  0.123456789;
    sctr->lon[1] =  sctr->lon[0] + 1.3235;
    sctr->lon[2] =  4.3445234;
    sctr->lon[3] =  sctr->lon[2] + 0.01;
    sctr->lon[4] =  0.0;
    sctr->lon[5] =  sctr->lon[4] + 2.0 * M_PI;

    /* Finally, we define spherical radii of the scattered cells.  Importantly,
     * the spherical radius is *constants* over the cells (but may vary from
     * cell to cell). */
    sctr->r[0]   =  shcs->r;
    sctr->r[1]   =  shcs->r + 1000.0;
    sctr->r[2]   =  shcs->r + 2000.0;
    /* --------------------------------------------------------------------- */


    /* Initialize an array to store the synthesized signal */
    f = (double *)malloc(grd->nlat * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis */
    charm_shs_cell(sctr, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Print the synthesized values */
    for (size_t i = 0; i < sctr->nlat; i++)
        printf("%0.16e\n", f[i]);


    charm_crd_free(sctr);
    free(f);


    printf("===========================\n\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Solid spherical harmonic synthesis at a grid "
           "of cells\n");
    printf("===========================\n");


    /* Initialize a "charm_crd" structure to hold a grid of cells with "15"
     * cells in the latitudinal direction and "30" cells in the longitudinal
     * direction. */
    grd = charm_crd_init(CHARM_CRD_CELLS_GRID, 15, 30);


    /* Define some grid cells */
    for (size_t i = 0; i < grd->nlat; i++)
    {
        /* The minimum cell latitudes */
        grd->lat[2 * i]      = -M_PI_2 + (double)i       * M_PI / grd->nlat;

        /* The maximum cell latitudes */
        grd->lat[2 * i + 1]  = -M_PI_2 + (double)(i + 1) * M_PI / grd->nlat;

        /* The spherical radii (may vary with the latitude, but we use
         * a constant value here, because of the example that follows after
         * this one). */
        grd->r[i] = shcs->r + 1000.0;
    }
    for (size_t j = 0; j < grd->nlon; j++)
    {
        /* The minimum cell longitudes */
        grd->lon[2 * j]      = (double)j       * (2.0 * M_PI) / grd->nlon;

        /* The maximum cell longitudes */
        grd->lon[2 * j + 1]  = (double)(j + 1) * (2.0 * M_PI) / grd->nlon;
    }


    /* Initialize an array to store the synthesized signal */
    f = (double *)malloc(grd->nlat * grd->nlon * sizeof(double));
    if (f == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the synthesis */
    charm_shs_cell(grd, shcs, nmax, f, err);
    charm_err_handler(err, 1);


    /* Print some of the synthesized values */
    /* --------------------------------------------------------------------- */
    i = 0;
    j = 10;
    printf("f(%zu, %zu) = %0.16e\n", i, j, f[i * grd->nlon + j]);


    i = 4;
    j = 3;
    printf("f(%zu, %zu) = %0.16e\n", i, j, f[i * grd->nlon + j]);
    /* --------------------------------------------------------------------- */


    printf("===========================\n\n");
    /* ===================================================================== */






    /* ===================================================================== */
    printf("\n");
    printf("Spherical harmonic analysis of block-mean values in cells\n");
    printf("===========================\n");


    /* In this example, we use the "grd" structure and the signal "f" from the
     * previous example. */


    /* Initialize a new structure of spherical harmonic coefficients for
     * coefficients to be recovered from the harmonic analysis.  Used are the
     * same scalling constants as in "shcs". */
    shcs2 = charm_shc_init(nmax, shcs->mu, shcs->r);
    if (shcs2 == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"charm_shc\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* Do the analysis using the method of the approximate quadrature */
    charm_sha_cell(grd, f, nmax, CHARM_SHA_CELL_AQ, shcs2, err);
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


    charm_crd_free(grd);
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
