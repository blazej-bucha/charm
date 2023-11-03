/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_cell_isSctr)(int type)
{
    return type == CHARM_CRD_CELL_SCATTERED;
}
