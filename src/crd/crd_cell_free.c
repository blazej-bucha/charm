/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(crd_cell_free)(CHARM(cell) *cell)
{
    if (cell == NULL)
        return;


    if (cell->owner)
    {
        free(cell->latmin);
        free(cell->latmax);
        free(cell->lonmin);
        free(cell->lonmax);
        free(cell->r);
    }
    free(cell);


    return;
}
