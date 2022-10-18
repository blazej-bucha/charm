/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */







int CHARM(crd_point_check_inputs)(int type, size_t nlat, size_t nlon)
{
    /* Check the inputs */
    /* --------------------------------------------------------------------- */
    /* Check "type" for supported values */
    if ((type != CHARM_CRD_POINTS_SCATTERED) &&
        (type != CHARM_CRD_POINTS_GRID) &&
        (type != CHARM_CRD_POINTS_GRID_GL) &&
        (type != CHARM_CRD_POINTS_GRID_DH1) &&
        (type != CHARM_CRD_POINTS_GRID_DH2))
        return 1;


    /* At least one point is required. */
    if (nlat < 1 || nlon < 1)
        return 2;


    /* For scattered points, "nlat" must be equal to "nlon" */
    if (type == CHARM_CRD_POINTS_SCATTERED)
    {
        if (nlat != nlon)
            return 3;
    }


    return 0;
    /* --------------------------------------------------------------------- */
}

