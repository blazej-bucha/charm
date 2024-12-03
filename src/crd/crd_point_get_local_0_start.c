/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_get_local_0_start.h"
/* ------------------------------------------------------------------------- */






size_t CHARM(crd_point_get_local_0_start)(const CHARM(point) *pnt)
{
#if HAVE_MPI
    return pnt->local_0_start;
#else
    return 0 * pnt->nlat;  /* This useless multiplication is just to suppress
                            * compiler warnings about not using "pnt" if
                            * "!HAVE_MPI". */
#endif
}
