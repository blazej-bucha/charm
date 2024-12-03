/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "shc_block_struct.h"
#include "shc_block_have_order.h"
/* ------------------------------------------------------------------------- */






/* Returns "1" if "shcs_block" contains coefficients of degree "m" and "0"
 * otherwise.
 *
 * This function can be called only after calling "charm_shc_block_init" and
 * "charm_shc_block_get_coeffs" (in that order). */
_Bool CHARM(shc_block_have_order)(const CHARM(shc_block) *shcs_block,
                                  unsigned long m)
{
    return ((m >= shcs_block->mfirst) && (m <= shcs_block->mlast));
}
