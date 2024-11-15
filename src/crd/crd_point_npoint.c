/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isSctr.h"
#include "crd_point_isGrid.h"
#include "crd_point_npoint.h"
/* ------------------------------------------------------------------------- */






size_t CHARM(crd_point_npoint)(int type,
                               size_t nlat,
                               size_t nlon)
{
    if (CHARM(crd_point_isSctr)(type))
        return nlat;
    else if (CHARM(crd_point_isGrid)(type))
        return nlat * nlon;
    else
        return 0;
}
