/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isGLGrid.h"
#include "crd_point_isDHGrid.h"
#include "crd_point_isQuadGrid.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_point_isQuadGrid)(int type)
{
    return CHARM(crd_point_isGLGrid)(type) || CHARM(crd_point_isDHGrid)(type);
}
