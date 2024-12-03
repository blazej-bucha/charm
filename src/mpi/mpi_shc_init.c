/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../shc/shc_init_chunk.h"
#include "../shc/shc_check_chunk_orders.h"
#include "mpi_count_aint.h"
#include "mpi_shc_check_struct.h"
#include "mpi_size_t2charm_mpi_count.h"
#include "mpi_allgather.h"
#include "mpi_gatherv.h"
#include "mpi_shc_local2distributed.h"
#include "mpi_err_gather.h"
#include "mpi_err_isdistributed.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(mpi_shc_init)(unsigned long nmax,
                                REAL mu,
                                REAL r,
                                REAL *c,
                                REAL *s,
                                size_t nchunk,
                                const unsigned long *order,
                                MPI_Comm comm,
                                CHARM(err) *err)
{
    /* ----------------------------------------------------------------- */
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    /* ----------------------------------------------------------------- */


    /* Memory allocations */
    /* ----------------------------------------------------------------- */
    CHARM(shc) *shcs            = NULL;
    CHARM_MPI_COUNT *nchunk_all = NULL;
    CHARM_MPI_COUNT *recv_cnt   = NULL;
    CHARM_MPI_AINT *displs      = NULL;
    unsigned long *order_all    = NULL;


    if (!CHARM(mpi_err_isdistributed)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }


    nchunk_all = (CHARM_MPI_COUNT *)malloc(size * sizeof(CHARM_MPI_COUNT));
    if (nchunk_all == NULL)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    recv_cnt = (CHARM_MPI_COUNT *)malloc(size * sizeof(CHARM_MPI_COUNT));
    if (recv_cnt == NULL)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    displs = (CHARM_MPI_AINT *)malloc(size * sizeof(CHARM_MPI_AINT));
    if (displs == NULL)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    if (!CHARM(mpi_err_isempty)(err))
        goto BARRIER;
    /* ----------------------------------------------------------------- */


    /* Prepare the output structure */
    /* ----------------------------------------------------------------- */
    shcs = CHARM(shc_init_chunk)(nmax, mu, r, c, s, nchunk, order, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }
    /* ----------------------------------------------------------------- */


    /* Check the distribution of chunks across MPI processes */
    /* ----------------------------------------------------------------- */
    /* Gather all "order" to a root */
    /* ................................................................. */
    /* Cast "nchunk" to "int" as required by MPI */
    CHARM_MPI_COUNT nchunk_int = CHARM(mpi_size_t2charm_mpi_count)(nchunk,
                                                                   err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }


BARRIER:
    /* Before the collective call, check that all previous calls where
     * successful */
    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;


    CHARM_MPI_ALLGATHER(&nchunk_int, 1, CHARM_MPI_COUNT_DATATYPE,
                        nchunk_all, 1, CHARM_MPI_COUNT_DATATYPE, comm);


    /* Sum the total number of chunks stored across all MPI processes */
    CHARM_MPI_COUNT nchunk_all_sum = 0;
    for (int i = 0; i < size; i++)
        nchunk_all_sum += nchunk_all[i];


    /* Array of the number of elements to be received from individual MPI
     * processes */
    for (int i = 0; i < size; i++)
        recv_cnt[i] = 2 * nchunk_all[i];  /* "times 2" because for each chunk
                                           * in "nchunk_all" we will receive
                                           * the minimum and the maximum orders
                                           * */


    /* Displacements to place "order" in "recv_cnt" */
    displs[0] = 0;
    for (int i = 1; i < size; i++)
        displs[i] = displs[i - 1] + 2 * (CHARM_MPI_AINT)nchunk_all[i - 1];


    /* "order_all" will store all values from "order" from all MPI processes,
     * with the exception being processes with "nchunk == 0" (that is, with no
     * chunks/data). */
    order_all = (unsigned long *)malloc(2 * nchunk_all_sum *
                                        sizeof(unsigned long));
    if (order_all == NULL)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    /* Another barrier before a collective call */
    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;


    CHARM_MPI_GATHERV(order, 2 * (int)nchunk, MPI_UNSIGNED_LONG,
                      order_all, recv_cnt, displs, MPI_UNSIGNED_LONG, 0, comm);
    /* ................................................................. */


    /* ................................................................. */
    if (rank == 0)
    {
        /* "shc_check_chunk_order" has already been called inside
         * "shc_init_chunk", but with "check_gaps" being false.  Now that we
         * have gathered boundaries of all chunks across MPI processes, we
         * shall call the function again, this time with "check_gaps = 1",
         * meaning that the function will now also check for possible gaps
         * between chunks. */
        CHARM(shc_check_chunk_orders)(nmax, (size_t)nchunk_all_sum, order_all,
                                      1, err);
        if (!CHARM(err_isempty)(err))
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            /* No "goto FAILURE" here.  If there was an error, all MPI
             * processes must go to "FAILURE", so we must wait until we are
             * outside the "if (rank == 0)". */
    }


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* ................................................................. */
    /* ----------------------------------------------------------------- */


    /* ----------------------------------------------------------------- */
    CHARM(mpi_shc_local2distributed)(shcs, comm);
    CHARM(mpi_shc_check_struct)(shcs, err);


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* ----------------------------------------------------------------- */


EXIT:
    free(nchunk_all);
    free(recv_cnt);
    free(displs);
    free(order_all);


    CHARM(mpi_err_gather)(err);


    return shcs;


FAILURE:
    CHARM(shc_free)(shcs);
    shcs = NULL;
    goto EXIT;
}
