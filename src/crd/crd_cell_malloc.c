/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_cell_alloc.h"
#include "crd_cell_check_inputs.h"
/* ------------------------------------------------------------------------- */







CHARM(cell) *CHARM(crd_cell_malloc)(int type, size_t nlat, size_t nlon)
{
    return CHARM(crd_cell_alloc)(type, nlat, nlon, malloc);
}

