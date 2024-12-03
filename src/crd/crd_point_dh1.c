/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_dh1_chunk.h"
#include "crd_point_quad.h"
/* ------------------------------------------------------------------------- */






CHARM(point) *CHARM(crd_point_dh1)(unsigned long nmax,
                                   REAL r)
{
    return CHARM(crd_point_quad)(nmax, r, CHARM(crd_point_dh1_shape),
                                 CHARM(crd_point_dh1_chunk));
}
