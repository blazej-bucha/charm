/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Increases the error level of an error structure.  This function is not
 * a part of the API. */
void CHARM(err_inc_level)(CHARM(err) *err)
{
    if ((err == NULL) || err->issaturated)
        return;


    /* Increase the level of the error structure */
    err->level += 1;


    /* If the "CHARM_ERR_MAX_LEVEL" has been reached, mark the error structure
     * as saturated to prevent further error propagation. */
    if (err->level >= CHARM_ERR_MAX_LEVEL)
        err->issaturated = 1;


    return;
}
