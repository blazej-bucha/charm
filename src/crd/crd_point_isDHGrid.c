/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isDH1Grid.h"
#include "crd_point_isDH2Grid.h"
#include "crd_point_isDHGrid.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_point_isDHGrid)(int type)
{
    return CHARM(crd_point_isDH1Grid)(type) ||
           CHARM(crd_point_isDH2Grid)(type);
}
