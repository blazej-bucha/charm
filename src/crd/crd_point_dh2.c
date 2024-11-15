/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_dh2_chunk.h"
#include "crd_point_quad.h"
/* ------------------------------------------------------------------------- */






CHARM(point) *CHARM(crd_point_dh2)(unsigned long nmax,
                                   REAL r)
{
    return CHARM(crd_point_quad)(nmax, r, CHARM(crd_point_dh2_shape),
                                 CHARM(crd_point_dh2_chunk));
}
