/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(crd_free)(CHARM(crd) *crd)
{
    if (crd == NULL)
        return;


    free(crd->lat);
    free(crd->lon);
    free(crd->r);
    free(crd->w);
    free(crd);


    return;
}
