/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isEmpty.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_point_isEmpty)(const CHARM(point) *pnt)
{
#if HAVE_MPI
    return pnt->local_npoint == 0;
#else
    return pnt->npoint == 0;
#endif
}
