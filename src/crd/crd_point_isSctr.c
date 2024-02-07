/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isSctr.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_point_isSctr)(int type)
{
    return type == CHARM_CRD_POINT_SCATTERED;
}
