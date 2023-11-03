/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_point_isGrid)(int type)
{
    return (type == CHARM_CRD_POINT_GRID) ||
           (type == CHARM_CRD_POINT_GRID_GL) ||
           (type == CHARM_CRD_POINT_GRID_DH1) ||
           (type == CHARM_CRD_POINT_GRID_DH2);
}
