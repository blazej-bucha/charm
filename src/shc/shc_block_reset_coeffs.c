/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <string.h>
#include "../prec.h"
#include "shc_block_struct.h"
#include "shc_block_reset_coeffs.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_block_reset_coeffs)(CHARM(shc_block) *shcs_block)
{
#if CHARM_OPENMP
#pragma omp master
#endif
    {
    memset(shcs_block->c, 0, shcs_block->ncs_max * sizeof(REAL));
    memset(shcs_block->s, 0, shcs_block->ncs_max * sizeof(REAL));
    }


    return;
}
