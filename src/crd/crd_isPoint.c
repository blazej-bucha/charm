/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isSctr.h"
#include "crd_point_isCustGrid.h"
#include "crd_point_isGLGrid.h"
#include "crd_point_isDHGrid.h"
#include "crd_isPoint.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_isPoint)(int type)
{
    return CHARM(crd_point_isSctr)(type) ||
           CHARM(crd_point_isCustGrid)(type) ||
           CHARM(crd_point_isGLGrid)(type) ||
           CHARM(crd_point_isDHGrid)(type);
}
