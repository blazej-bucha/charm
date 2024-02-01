/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_cell_isGrid.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_cell_isGrid)(int type)
{
    return type == CHARM_CRD_CELL_GRID;
}
