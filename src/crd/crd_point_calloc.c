/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_alloc.h"
#include "../misc/misc_calloc.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(crd_point_calloc)(int type, size_t nlat, size_t nlon)
{
    return CHARM(crd_point_alloc)(type, nlat, nlon, CHARM(misc_calloc));
}

