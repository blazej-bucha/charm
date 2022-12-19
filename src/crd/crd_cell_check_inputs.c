/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */







int CHARM(crd_cell_check_inputs)(int type, size_t nlat, size_t nlon)
{
    /* Check the inputs */
    /* --------------------------------------------------------------------- */
    /* Check "type" for supported values */
    if ((type != CHARM_CRD_CELL_SCATTERED) && (type != CHARM_CRD_CELL_GRID))
        return 1;


    /* At least one cell is required. */
    if (nlat < 1 || nlon < 1)
        return 2;


    /* For scattered cells, "nlat" must be equal to "nlon" */
    if (type == CHARM_CRD_CELL_SCATTERED)
    {
        if (nlat != nlon)
            return 3;
    }


    return 0;
    /* --------------------------------------------------------------------- */
}

