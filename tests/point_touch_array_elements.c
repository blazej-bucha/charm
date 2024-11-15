/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/crd/crd_point_isQuadGrid.h"
#include "../src/crd/crd_point_get_local_nlat.h"
#include "../src/crd/crd_point_get_local_nlon.h"
#include "point_touch_array_elements.h"
/* ------------------------------------------------------------------------- */






/* Check that "pnt->lat", "pnt->lon", "pnt->r" and possibly also "pnt->w" have
 * access to the required number of elements.  The bad thing is that if this
 * condition is not satisfied, a segfault error will appear on runtime, not
 * leaving any room to report to the user what kind of error was encountered.
 * Still, in this rare case, a segfault error is still much better than no
 * error at all. */
/* ------------------------------------------------------------------------- */
void point_touch_array_elements(CHARM(point) *pnt)
{
    if (pnt == NULL)
        return;


    size_t local_nlat = CHARM(crd_point_get_local_nlat)(pnt);
    size_t local_nlon = CHARM(crd_point_get_local_nlon)(pnt);


    for (size_t i = 0; i < local_nlat; i++)
        pnt->lat[i] = PREC(1.0);
    for (size_t i = 0; i < local_nlon; i++)
        pnt->lon[i] = PREC(1.0);
    for (size_t i = 0; i < local_nlat; i++)
        pnt->r[i] = PREC(1.0);
    if (CHARM(crd_point_isQuadGrid)(pnt->type))
    {
        for (size_t i = 0; i < local_nlat; i++)
            pnt->w[i] = PREC(1.0);
    }


    return;
}
/* ------------------------------------------------------------------------- */
