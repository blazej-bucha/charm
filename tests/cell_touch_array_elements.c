/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cell_touch_array_elements.h"
/* ------------------------------------------------------------------------- */






/* Check that "cell->latmin", etc. have access to the required number of
 * elements.  The bad thing is that if this condition is not satisfied,
 * a segfault error will appear on runtime, not leaving any room to report to
 * the user what kind of error was encountered.  Still, in this rare case,
 * a segfault error is still much better than no error at all. */
/* ------------------------------------------------------------------------- */
void cell_touch_array_elements(CHARM(cell) *cell)
{
    if (cell == NULL)
        return;


    for (size_t i = 0; i < cell->nlat; i++)
    {
        cell->latmin[i] = PREC(1.0);
        cell->latmax[i] = PREC(1.0);
    }
    for (size_t i = 0; i < cell->nlon; i++)
    {
        cell->lonmin[i] = PREC(1.0);
        cell->lonmax[i] = PREC(1.0);
    }
    for (size_t i = 0; i < cell->nlat; i++)
        cell->r[i] = PREC(1.0);


    return;
}
/* ------------------------------------------------------------------------- */
