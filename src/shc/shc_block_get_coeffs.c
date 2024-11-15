/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <string.h>
#include "../prec.h"
#if HAVE_MPI
#   include <mpi.h>
#   include "../mpi/mpi_real.h"
#   include "../mpi/mpi_count_aint.h"
#   include "../mpi/mpi_ibcast.h"
#endif
#include "../err/err_propagate.h"
#include "../err/err_isempty_all_mpi_processes.h"
#include "shc_block_struct.h"
#include "shc_block_set_mfirst.h"
#include "shc_block_get_coeffs.h"
/* ------------------------------------------------------------------------- */






/* Copies a block of coefficients of orders "mfirst", "mfirst + 1", ... from
 * the process, which stores these coefficients in "shcs", to "shcs_block" on
 * all processes.  The total number of orders that are gathered depends mostly
 * on the maximum order of the chunk, where "mfirst" is found. */
void CHARM(shc_block_get_coeffs)(const CHARM(shc) *shcs
#if HAVE_MPI
                                 , CHARM(shc_block) *shcs_block,
                                 unsigned long mfirst,
                                 CHARM(err) *err
#endif
                                )
{
    /* There is nothing we can do for non-distributed "shcs" */
    /* --------------------------------------------------------------------- */
    if (!shcs->distributed)
        return;
    /* --------------------------------------------------------------------- */


#if HAVE_MPI
    /* This code is for distributed "charm_shc" structures only */
#   if HAVE_OPENMP
#pragma omp master
#   endif
    {
    /* Set "mfirst" of "shcs_block" and update the rest of the metadata of
     * "shcs_block" */
    /* ----------------------------------------------------------------- */
    CHARM(shc_block_set_mfirst)(shcs_block, shcs, mfirst, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto EXIT_STRUCTURED_BLOCK;
    /* ----------------------------------------------------------------- */


    /* Make a hard copy of the coefficients on the "shcs_block->root" */
    /* ----------------------------------------------------------------- */
    int rank;
    MPI_Comm_rank(shcs_block->comm, &rank);


    if (rank == shcs_block->root)
    {
        memcpy(shcs_block->c, shcs->c[mfirst], shcs_block->nc * sizeof(REAL));
        memcpy(shcs_block->s, shcs->s[mfirst], shcs_block->ns * sizeof(REAL));
    }
    /* ----------------------------------------------------------------- */


    /* Now broadcast the coefficients from the "shc_block->root" */
    /* ----------------------------------------------------------------- */
    MPI_Request reqs[2];
    CHARM_MPI_IBCAST(shcs_block->c, (CHARM_MPI_COUNT)shcs_block->nc,
                     CHARM_MPI_REAL, shcs_block->root, shcs_block->comm,
                     &reqs[0]);
    CHARM_MPI_IBCAST(shcs_block->s, (CHARM_MPI_COUNT)shcs_block->ns,
                     CHARM_MPI_REAL, shcs_block->root, shcs_block->comm,
                     &reqs[1]);
    MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);
    /* ----------------------------------------------------------------- */


EXIT_STRUCTURED_BLOCK:
    ;
    }
#   if HAVE_OPENMP
#pragma omp barrier
#   endif
#endif


    return;
}
