/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_cell_check_inputs.h"
/* ------------------------------------------------------------------------- */







CHARM(cell) *CHARM(crd_cell_init)(int type, size_t nlat, size_t nlon,
                                  REAL *latmin, REAL *latmax,
                                  REAL *lonmin, REAL *lonmax,
                                  REAL *r)
{
    /* Check the inputs */
    /* --------------------------------------------------------------------- */
    if (CHARM(crd_cell_check_inputs)(type, nlat, nlon))
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Initialize the "CHARM(cell)" structure */
    /* --------------------------------------------------------------------- */
    /* Allocate memory for the "CHARM(cell)" data type */
    CHARM(cell) *cell = (CHARM(cell) *)malloc(sizeof(CHARM(cell)));
    if (cell == NULL)
        return cell;
    /* --------------------------------------------------------------------- */


    /* Set the array members of "cell" */
    /* --------------------------------------------------------------------- */
    cell->latmin = latmin;
    cell->latmax = latmax;
    cell->lonmin = lonmin;
    cell->lonmax = lonmax;
    cell->r      = r;
    /* --------------------------------------------------------------------- */


    /* Set the scalar members of "cell" */
    /* --------------------------------------------------------------------- */
    cell->nlat  = nlat;
    cell->nlon  = nlon;
    if (type == CHARM_CRD_CELL_SCATTERED)
        cell->ncell = nlat;
    else
        cell->ncell = nlat * nlon;
    cell->type  = type;
    cell->owner = 0;
    /* --------------------------------------------------------------------- */


    return cell;
}
