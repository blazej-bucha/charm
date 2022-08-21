/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <string.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* An internal function to set all coefficients of a "shc" structure to
 * zero. */
void CHARM(shc_reset_coeffs)(CHARM(shc) *shcs)
{
    memset(shcs->c[0], 0, shcs->nc * sizeof(REAL));
    memset(shcs->s[0], 0, shcs->ns * sizeof(REAL));


    return;
}
