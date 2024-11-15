/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_point_check_inputs.h"
#include "crd_point_isQuadGrid.h"
#include "crd_point_alloc.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(crd_point_alloc)(int type,
                                     size_t nlat,
                                     size_t nlon,
                                     void *(*alloc)(size_t))
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


    if (nlat > 0)
    {
        lat = (REAL *)alloc(nlat * sizeof(REAL));
        if (lat == NULL)
            goto FAILURE;


        r = (REAL *)alloc(nlat * sizeof(REAL));
        if (r == NULL)
            goto FAILURE;
    }


    if (nlon > 0)
    {
        lon = (REAL *)alloc(nlon * sizeof(REAL));
        if (lon == NULL)
            goto FAILURE;
    }


    pnt = CHARM(crd_point_init)(type, nlat, nlon, lat, lon, r);
    if (pnt == NULL)
        goto FAILURE;


    if (CHARM(crd_point_isQuadGrid)(type))
    {
        w = (REAL *)alloc(nlat * sizeof(REAL));
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

