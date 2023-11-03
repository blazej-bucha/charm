/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../crd/crd_isCell.h"
#include "../crd/crd_cell_isSctr.h"
/* ------------------------------------------------------------------------- */







int CHARM(crd_cell_check_inputs)(int type, size_t nlat, size_t nlon)
{
    /* Check the inputs */
    /* --------------------------------------------------------------------- */
    /* Check "type" for supported values */
    if (!CHARM(crd_isCell)(type))
        return 1;


    /* At least one cell is required. */
    if (nlat < 1 || nlon < 1)
        return 2;


    /* For scattered cells, "nlat" must be equal to "nlon" */
    if (CHARM(crd_cell_isSctr)(type))
    {
        if (nlat != nlon)
            return 3;
    }


    return 0;
    /* --------------------------------------------------------------------- */
}

