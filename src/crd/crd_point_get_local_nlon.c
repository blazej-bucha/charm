/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_get_local_nlon.h"
/* ------------------------------------------------------------------------- */






/* Returns the number of longitudes, for which the synthesis/analysis with
 * point data is to be performed locally. */
size_t CHARM(crd_point_get_local_nlon)(const CHARM(point) *pnt)
{
#if HAVE_MPI
    return pnt->local_nlon;
#else
    return pnt->nlon;
#endif
}
