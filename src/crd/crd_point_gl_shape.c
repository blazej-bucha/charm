/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_quad_l.h"
/* ------------------------------------------------------------------------- */






void CHARM(crd_point_gl_shape)(unsigned long nmax,
                               size_t *nlat,
                               size_t *nlon)
{
    unsigned long L = CHARM(crd_point_quad_l)(nmax);
    *nlat = L;
    *nlon = 2 * L;


    return;
}
