/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../misc/misc_nan.h"
#include "shc_block_struct.h"
#include "shc_block_nan.h"
/* ------------------------------------------------------------------------- */






/* Sets all the coefficients of "shcs_block" to "NAN" to indicate an error.  As
 * a consequence, the "NAN"s will be propagated to the final result, thereby
 * indicating the error. */
void CHARM(shc_block_nan)(CHARM(shc_block) *shcs_block)
{
    for (size_t i = 0; i < shcs_block->ncs_max; i++)
        shcs_block->c[i] = NAN;
    for (size_t i = 0; i < shcs_block->ncs_max; i++)
        shcs_block->s[i] = NAN;


    return;
}
