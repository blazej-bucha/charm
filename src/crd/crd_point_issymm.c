/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../crd/crd_point_isEmpty.h"
#include "../crd/crd_point_isGrid.h"
#include "../misc/misc_is_nearly_equal.h"
#include "crd_point_issymm.h"
/* ------------------------------------------------------------------------- */






_Bool CHARM(crd_point_issymm)(const CHARM(point) *pnt)
{
    if (CHARM(crd_point_isEmpty)(pnt))
        return 0;


    if (!CHARM(crd_point_isGrid)(pnt->type))
        return 0;


    size_t local_nlat;
#if HAVE_MPI
    local_nlat = pnt->local_nlat;
#else
    local_nlat = pnt->nlat;
#endif
    if (local_nlat < 2)
        return 0;


    /* Now "pnt" has to be a non-empty grid with more than 1 latitude */


    /* Check locally that each positive latitude has its negative counterpart.
     * We do this also with the zero latitude if present. */
    for (size_t i = 0; i < (local_nlat + 1) / 2; i++)
        if (!CHARM(misc_is_nearly_equal)(FABS(pnt->lat[i]),
                                         FABS(pnt->lat[local_nlat - 1 - i]),
                                         CHARM(glob_threshold)))
            return 0;


    /* If the local number of latitudes is odd and is larger than "1", the
     * latitude in the middle must be zero for a grid to be symmetric */
    if (local_nlat % 2)
        if (!CHARM(misc_is_nearly_equal)(pnt->lat[local_nlat / 2],
                                         PREC(0.0),
                                         CHARM(glob_threshold)))
            return 0;


    return 1;
}
