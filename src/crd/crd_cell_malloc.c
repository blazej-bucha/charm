/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_cell_check_inputs.h"
/* ------------------------------------------------------------------------- */







CHARM(cell) *CHARM(crd_cell_malloc)(int type, size_t nlat, size_t nlon)
{
    /* Check the inputs */
    /* --------------------------------------------------------------------- */
    if (CHARM(crd_cell_check_inputs)(type, nlat, nlon))
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Set the members of "cell" */
    /* --------------------------------------------------------------------- */
    CHARM(cell) *cell = NULL;
    REAL *latmin      = NULL;
    REAL *latmax      = NULL;
    REAL *lonmin      = NULL;
    REAL *lonmax      = NULL;
    REAL *r           = NULL;


    latmin = (REAL *)malloc(nlat * sizeof(REAL));
    if (latmin == NULL)
        goto FAILURE;


    latmax = (REAL *)malloc(nlat * sizeof(REAL));
    if (latmax == NULL)
        goto FAILURE;


    lonmin = (REAL *)malloc(nlon * sizeof(REAL));
    if (lonmin == NULL)
        goto FAILURE;


    lonmax = (REAL *)malloc(nlon * sizeof(REAL));
    if (lonmax == NULL)
        goto FAILURE;


    r = (REAL *)malloc(nlat * sizeof(REAL));
    if (r == NULL)
        goto FAILURE;


    cell = CHARM(crd_cell_init)(type, nlat, nlon, latmin, latmax,
                                lonmin, lonmax, r);
    if (cell == NULL)
        goto FAILURE;


    cell->owner = 1;
    /* --------------------------------------------------------------------- */


EXIT:
    /* --------------------------------------------------------------------- */
    return cell;
    /* --------------------------------------------------------------------- */


FAILURE:
    /* --------------------------------------------------------------------- */
    free(latmin);
    free(latmax);
    free(lonmin);
    free(lonmax);
    free(r);
    free(cell);


    cell = NULL;


    goto EXIT;
    /* --------------------------------------------------------------------- */
}

