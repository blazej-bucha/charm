/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_free)(CHARM(shc) *shcs)
{
    if (shcs == NULL)
        return;


    free(shcs->c[0]);
    free(shcs->s[0]);
    free(shcs->c);
    free(shcs->s);
    free(shcs);


    return;
}
