/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_isPoint.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_isPoint)(int type)
{
    return (type == CHARM_CRD_POINT_SCATTERED) ||
           (type == CHARM_CRD_POINT_GRID) ||
           (type == CHARM_CRD_POINT_GRID_GL) ||
           (type == CHARM_CRD_POINT_GRID_DH1) ||
           (type == CHARM_CRD_POINT_GRID_DH2);
}
