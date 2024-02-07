/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_isCell.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_isCell)(int type)
{
    return (type == CHARM_CRD_CELL_SCATTERED) ||
           (type == CHARM_CRD_CELL_GRID);
}
