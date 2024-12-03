/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"
#include "crd_point_check_inputs.h"
#include "crd_point_isSctr.h"
#include "crd_point_isGrid.h"
#include "crd_point_npoint.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(crd_point_init)(int type,
                                    size_t nlat,
                                    size_t nlon,
                                    REAL *lat,
                                    REAL *lon,
                                    REAL *r)
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
    pnt->lat = pnt->lon = pnt->r = pnt->w = NULL;


    if (nlat > 0)
    {
        pnt->lat = lat;
        pnt->r   = r;
    }


    if (nlon > 0)
        pnt->lon = lon;
    /* --------------------------------------------------------------------- */


    /* Set the scalar members of "pnt" */
    /* --------------------------------------------------------------------- */
    pnt->nlat = nlat;
    pnt->nlon = nlon;
    pnt->npoint = CHARM(crd_point_npoint)(type, nlat, nlon);
    pnt->type = type;
    pnt->owner = 0;
    pnt->distributed = 0;
#if HAVE_MPI
    pnt->local_nlat = pnt->nlat;
    pnt->local_nlon = pnt->nlon;
    pnt->local_npoint = pnt->npoint;
    pnt->local_0_start = 0;
    pnt->comm = MPI_COMM_NULL;
#endif
    /* --------------------------------------------------------------------- */


    return pnt;
}
