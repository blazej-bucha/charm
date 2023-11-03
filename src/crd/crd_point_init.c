/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_point_check_inputs.h"
#include "crd_point_isSctr.h"
#include "crd_point_isGrid.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(crd_point_init)(int type, size_t nlat, size_t nlon,
                                    REAL *lat, REAL *lon, REAL *r)
{
    /* Check the inputs */
    /* --------------------------------------------------------------------- */
    if (CHARM(crd_point_check_inputs)(type, nlat, nlon))
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Initialize the "CHARM(point)" structure */
    /* --------------------------------------------------------------------- */
    /* Allocate memory for the "CHARM(point)" data type */
    CHARM(point) *pnt = (CHARM(point) *)malloc(sizeof(CHARM(point)));
    if (pnt == NULL)
        return pnt;
    /* --------------------------------------------------------------------- */


    /* Set the array members of "pnt" */
    /* --------------------------------------------------------------------- */
    pnt->lat = lat;
    pnt->lon = lon;
    pnt->r   = r;
    pnt->w   = NULL;
    /* --------------------------------------------------------------------- */


    /* Set the scalar members of "pnt" */
    /* --------------------------------------------------------------------- */
    pnt->nlat  = nlat;
    pnt->nlon  = nlon;
    if (CHARM(crd_point_isSctr)(type))
        pnt->npoint = nlat;
    else if (CHARM(crd_point_isGrid)(type))
        pnt->npoint = nlat * nlon;
    pnt->type  = type;
    pnt->owner = 0;
    /* --------------------------------------------------------------------- */


    return pnt;
}
