/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_isPoint.h"
#include "crd_point_isSctr.h"
#include "crd_point_check_inputs.h"
/* ------------------------------------------------------------------------- */







int CHARM(crd_point_check_inputs)(int type,
                                  size_t nlat,
                                  size_t nlon)
{
    /* Check "type" for supported values */
    if (!CHARM(crd_isPoint)(type))
        return 1;


    /* For scattered points, "nlat" must be equal to "nlon" */
    if (CHARM(crd_point_isSctr)(type))
    {
        if (nlat != nlon)
            return 2;
    }


    return 0;
}

