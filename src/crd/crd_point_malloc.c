/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_point_alloc.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(crd_point_malloc)(int type, size_t nlat, size_t nlon)
{
    return CHARM(crd_point_alloc)(type, nlat, nlon, malloc);
}

