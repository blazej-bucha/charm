/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */







CHARM(crd) *CHARM(crd_init)(int type, size_t nlat, size_t nlon)
{
    /* Check the inputs */
    /* --------------------------------------------------------------------- */
    /* Check "type" for supported values */
    if ((type != CHARM_CRD_CELLS_GRID) &&
        (type != CHARM_CRD_CELLS_SCATTERED) &&
        (type != CHARM_CRD_POINTS_SCATTERED) &&
        (type != CHARM_CRD_POINTS_GRID) &&
        (type != CHARM_CRD_POINTS_GRID_GL) &&
        (type != CHARM_CRD_POINTS_GRID_DH1) &&
        (type != CHARM_CRD_POINTS_GRID_DH2))
        return NULL;


    /* At least one cell/point is required. */
    if (nlat < 1 || nlon < 1)
        return NULL;


    /* For scattered cells/points, "nlat" must be equal to "nlon" */
    if ((type == CHARM_CRD_CELLS_SCATTERED) ||
        (type == CHARM_CRD_POINTS_SCATTERED))
    {
        if (nlat != nlon)
            return NULL;
    }
    /* --------------------------------------------------------------------- */


    /* Initialize the "CHARM(crd)" structure */
    /* --------------------------------------------------------------------- */
    /* Allocate memory for the "CHARM(crd)" data type */
    CHARM(crd) *crd = (CHARM(crd) *)malloc(sizeof(CHARM(crd)));
    if (crd == NULL)
        return NULL;


    crd->lat = crd->lon = crd->w = crd->r = NULL;
    /* --------------------------------------------------------------------- */


    /* Allocate and initialize latitudes and longitudes */
    /* --------------------------------------------------------------------- */
    /* At first, get the number of latitudes and longitudes to be initialized
     * (depend on "type") */
    size_t nlat_all, nlon_all;
    if ((type == CHARM_CRD_CELLS_SCATTERED) ||
        (type == CHARM_CRD_CELLS_GRID))
    {
        /* Cells */
        nlat_all = 2 * nlat;
        nlon_all = 2 * nlon;
    }
    else if ((type != CHARM_CRD_POINTS_SCATTERED) ||
             (type != CHARM_CRD_POINTS_GRID) ||
             (type != CHARM_CRD_POINTS_GRID_GL) ||
             (type != CHARM_CRD_POINTS_GRID_DH1) ||
             (type != CHARM_CRD_POINTS_GRID_DH2))
    {
        /* Points */
        nlat_all = nlat;
        nlon_all = nlon;
    }


    /* Now do the allocation and initializations */
    crd->lat = (REAL *)calloc(nlat_all, sizeof(REAL));
    if (crd->lat == NULL)
        goto FAILURE;


    crd->lon = (REAL *)calloc(nlon_all, sizeof(REAL));
    if (crd->lon == NULL)
        goto FAILURE;


    crd->r = (REAL *)calloc(nlat, sizeof(REAL));
    if (crd->r == NULL)
        goto FAILURE;


    /* The following point grids are associated also with integration weights,
     * so let's allocate and initialize a "w" array to store the weights. */
    if ((type == CHARM_CRD_POINTS_GRID_GL) ||
        (type == CHARM_CRD_POINTS_GRID_DH1) ||
        (type == CHARM_CRD_POINTS_GRID_DH2))
    {
        crd->w = (REAL *)calloc(nlat, sizeof(REAL));
        if (crd->w == NULL)
            goto FAILURE;
    }
    /* --------------------------------------------------------------------- */


    /* Assign the number of latitudes and longitudes and the "type" value to
     * the "crd" struct */
    /* --------------------------------------------------------------------- */
    crd->nlat = nlat;
    crd->nlon = nlon;
    crd->type = type;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
EXIT:
    return crd;
    /* --------------------------------------------------------------------- */


    /* Release internally allocated memory if any memory allocation failed. */
    /* --------------------------------------------------------------------- */
FAILURE:
    free(crd->lat);
    free(crd->lon);
    free(crd->r);
    free(crd->w);
    free(crd);


    crd = NULL;


    goto EXIT;
    /* --------------------------------------------------------------------- */
}
