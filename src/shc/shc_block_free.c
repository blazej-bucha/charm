/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shc_block_struct.h"
#include "shc_block_free.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_block_free)(CHARM(shc_block) *shcs_block)
{
    if (shcs_block == NULL)
        return;


    if (shcs_block->owner)
    {
        free(shcs_block->c);
        free(shcs_block->s);
    }
#if HAVE_MPI
    free(shcs_block->have_m_all);  /* Always freed, regardless of
                                    * "shcs_block->owner" */
#endif
    free(shcs_block);


    return;
}
