/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "shc_block_struct.h"
#include "shc_block_get_idx.h"
/* ------------------------------------------------------------------------- */






/* Returns an index of "C_{mcurr, mcurr}" and "S_{mcurr, mcurr}" in
 * "shcs_block->c" and "shcs_block->c" with respect to the chunk starting with
 * order "shcs_block_mfirst" with "mcurr >= shcs_block_mfirst". */
unsigned long CHARM(shc_block_get_idx)(const CHARM(shc_block) *shcs_block,
                                       unsigned long mcurr)
{
    unsigned long nmaxp2 = shcs_block->nmax + 2;
    unsigned long idx = 0;
    for (unsigned long m = shcs_block->mfirst + 1; m <= mcurr; m++)
        idx += nmaxp2 - m;


    return idx;
}
