/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_propagate.h"
#include "../err/err_isempty_all_mpi_processes.h"
#include "shc_block_struct.h"
#include "shc_block_get_mlast_ncs_root.h"
#include "shc_block_set_mfirst.h"
/* ------------------------------------------------------------------------- */






/* Sets "shcs_block->mfirst" to "mfirst" and updates the rest of the metadata
 * of "shcs_block" */
void CHARM(shc_block_set_mfirst)(CHARM(shc_block) *shcs_block,
                                 const CHARM(shc) *shcs,
                                 unsigned long mfirst,
                                 CHARM(err) *err)
{
#if CHARM_OPENMP
#pragma omp master
#endif
    {
    unsigned long mlast;
    size_t ncs;
    int root;
    CHARM(shc_block_get_mlast_ncs_root)(shcs, shcs_block, mfirst, &mlast, &ncs,
                                        &root, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto EXIT_STRUCTURED_BLOCK;


    shcs_block->mfirst = mfirst;
    shcs_block->mlast  = mlast;
    shcs_block->nc     = ncs;
    shcs_block->ns     = ncs;
    shcs_block->root   = root;


EXIT_STRUCTURED_BLOCK:
    ;
    }


    return;
}
