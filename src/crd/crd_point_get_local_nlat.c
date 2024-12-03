/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_get_local_nlat.h"
/* ------------------------------------------------------------------------- */






/* Returns the number of latitudes, for which the synthesis/analysis with point
 * data is to be performed locally. */
size_t CHARM(crd_point_get_local_nlat)(const CHARM(point) *pnt)
{
#if HAVE_MPI
    return pnt->local_nlat;
#else
    return pnt->nlat;
#endif
}
