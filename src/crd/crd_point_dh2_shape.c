/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_quad_l.h"
/* ------------------------------------------------------------------------- */






void CHARM(crd_point_dh2_shape)(unsigned long nmax,
                                size_t *nlat,
                                size_t *nlon)
{
    unsigned long L = CHARM(crd_point_quad_l)(nmax);
    *nlat = 2 * L;
    *nlon = 4 * L;


    return;
}
