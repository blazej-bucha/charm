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
#   include "../mpi/mpi_ireduce.h"
#endif
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_isempty_all_mpi_processes.h"
#include "shc_block_struct.h"
#include "shc_block_get_mlast_ncs_root.h"
#include "shc_block_reset_coeffs.h"
#include "shc_block_set_coeffs.h"
/* ------------------------------------------------------------------------- */






/* Sums coefficients from "shcs_block" across all MPI processes and stores the
 * results at the MPI processes holding the chunk of coefficients given by
 * "mfirst" and "mlast". */
void CHARM(shc_block_set_coeffs)(CHARM(shc) *shcs
#if HAVE_MPI
                                 , CHARM(shc_block) *shcs_block,
                                 unsigned long mfirst,
                                 unsigned long mlast,
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
#   if CHARM_OPENMP
#pragma omp master
#   endif
    {
    /* --------------------------------------------------------------------- */
    REAL *c = NULL;
    REAL *s = NULL;


    c = (REAL *)malloc(shcs_block->nc * sizeof(REAL));
    if (c == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto BARRIER;
    }


    s = (REAL *)malloc(shcs_block->ns * sizeof(REAL));
    if (s == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto BARRIER;
    }


BARRIER:
    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto EXIT_STRUCTURED_BLOCK;


    MPI_Request reqs[2];
    CHARM_MPI_IREDUCE(shcs_block->c, c, (CHARM_MPI_COUNT)shcs_block->nc,
                      CHARM_MPI_REAL, MPI_SUM, shcs_block->root,
                      shcs_block->comm, &reqs[0]);
    CHARM_MPI_IREDUCE(shcs_block->s, s, (CHARM_MPI_COUNT)shcs_block->ns,
                     CHARM_MPI_REAL, MPI_SUM, shcs_block->root,
                     shcs_block->comm, &reqs[1]);
    MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);


    int rank;
    MPI_Comm_rank(shcs_block->comm, &rank);
    if (rank == shcs_block->root)
    {
        for (size_t i = 0; i < shcs_block->nc; i++)
            shcs->c[mfirst][i] += c[i];


        for (size_t i = 0; i < shcs_block->ns; i++)
            shcs->s[mfirst][i] += s[i];
    }
    /* --------------------------------------------------------------------- */


    /* Update some data in "shcs_block", so that the structure is ready for use
     * in the next iteration.  This needs to be done, however, only if the
     * maximum degree of "shcs" will not be exceeded. */
    /* --------------------------------------------------------------------- */
    unsigned long mfirst_new = mlast + 1;


    if (mfirst_new <= shcs->nmax)
    {
        unsigned long mlast_new;
        int root;
        size_t ncs;
        CHARM(shc_block_get_mlast_ncs_root)(shcs, shcs_block, mfirst_new,
                                            &mlast_new, &ncs, &root, err);
        if (!CHARM(err_isempty)(err))
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


        if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
            goto EXIT_STRUCTURED_BLOCK;


        shcs_block->mfirst = mfirst_new;
        shcs_block->mlast  = mlast_new;
        shcs_block->nc     = ncs;
        shcs_block->ns     = ncs;
        shcs_block->root   = root;
    }
    /* --------------------------------------------------------------------- */


    /* Reset the coefficients in "shcs_block" to zero, so that "shcs_block" can
     * be used with the next iteration */
    CHARM(shc_block_reset_coeffs)(shcs_block);


EXIT_STRUCTURED_BLOCK:
    free(c);
    free(s);
    }
#   if CHARM_OPENMP
#pragma omp barrier
#   endif
#endif
    /* --------------------------------------------------------------------- */


    return;
}
