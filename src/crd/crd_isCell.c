/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_cell_isGrid.h"
#include "crd_cell_isSctr.h"
#include "crd_isCell.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_isCell)(int type)
{
    return CHARM(crd_cell_isGrid)(type) || CHARM(crd_cell_isSctr)(type);
}
