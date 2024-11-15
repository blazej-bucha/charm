/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_gl_chunk.h"
#include "crd_point_quad.h"
/* ------------------------------------------------------------------------- */






CHARM(point) *CHARM(crd_point_gl)(unsigned long nmax,
                                  REAL r)
{
    return CHARM(crd_point_quad)(nmax, r, CHARM(crd_point_gl_shape),
                                 CHARM(crd_point_gl_chunk));
}
