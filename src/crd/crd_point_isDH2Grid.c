/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isDH2Grid.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_point_isDH2Grid)(int type)
{
    return type == CHARM_CRD_POINT_GRID_DH2;
}
