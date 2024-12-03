/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isGLGrid.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_point_isGLGrid)(int type)
{
    return type == CHARM_CRD_POINT_GRID_GL;
}
