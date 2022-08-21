/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(crd_point_free)(CHARM(point) *pnt)
{
    if (pnt == NULL)
        return;


    if (pnt->owner)
    {
        free(pnt->lat);
        free(pnt->lon);
        free(pnt->r);
        free(pnt->w);
    }
    free(pnt);


    return;
}
