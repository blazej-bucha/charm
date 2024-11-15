/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_isCell.h"
#include "crd_cell_isSctr.h"
#include "crd_cell_check_inputs.h"
/* ------------------------------------------------------------------------- */







int CHARM(crd_cell_check_inputs)(int type,
                                 size_t nlat,
                                 size_t nlon)
{
    /* Check the inputs */
    /* --------------------------------------------------------------------- */
    /* Check "type" for supported values */
    if (!CHARM(crd_isCell)(type))
        return 1;


    /* For scattered cells, "nlat" must be equal to "nlon" */
    if (CHARM(crd_cell_isSctr)(type))
    {
        if (nlat != nlon)
            return 2;
    }


    return 0;
    /* --------------------------------------------------------------------- */
}

