/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isCustGrid.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_point_isCustGrid)(int type)
{
    return type == CHARM_CRD_POINT_GRID;
}
