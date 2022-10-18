/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_point_check_inputs.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(crd_point_calloc)(int type, size_t nlat, size_t nlon)
{
    /* Check the inputs */
    /* --------------------------------------------------------------------- */
    if (CHARM(crd_point_check_inputs)(type, nlat, nlon))
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Set the members of "pnt" */
    /* --------------------------------------------------------------------- */
    CHARM(point) *pnt = NULL;
    REAL *lat         = NULL;
    REAL *lon         = NULL;
    REAL *r           = NULL;
    REAL *w           = NULL;


    lat = (REAL *)calloc(nlat, sizeof(REAL));
    if (lat == NULL)
        goto FAILURE;


    lon = (REAL *)calloc(nlon, sizeof(REAL));
    if (lon == NULL)
        goto FAILURE;


    r = (REAL *)calloc(nlat, sizeof(REAL));
    if (r == NULL)
        goto FAILURE;


    pnt = CHARM(crd_point_init)(type, nlat, nlon, lat, lon, r);
    if (pnt == NULL)
        goto FAILURE;


    if ((type == CHARM_CRD_POINTS_GRID_GL) ||
        (type == CHARM_CRD_POINTS_GRID_DH1) ||
        (type == CHARM_CRD_POINTS_GRID_DH2))
    {
        w = (REAL *)calloc(nlat, sizeof(REAL));
        if (w == NULL)
            goto FAILURE;

        pnt->w = w;
    }


    pnt->owner = 1;
    /* --------------------------------------------------------------------- */


EXIT:
    /* --------------------------------------------------------------------- */
    return pnt;
    /* --------------------------------------------------------------------- */


FAILURE:
    /* --------------------------------------------------------------------- */
    free(lat);
    free(lon);
    free(r);
    free(w);
    free(pnt);


    pnt = NULL;


    goto EXIT;
    /* --------------------------------------------------------------------- */
}

