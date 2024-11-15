/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../src/prec.h"
#include "error_messages.h"
#include "partition_interval.h"
#include "mpi_shc_distribute.h"
/* ------------------------------------------------------------------------- */






/* This function is not completely general.  To make it work, "nmax" must be
 * larger enough than "size". */
unsigned long *mpi_shc_distribute(unsigned long nmax,
                                  int rank,
                                  int size,
                                  size_t local_nchunk)
{
    /* Now split "[start, stop]" evenly across "local_nchunk" chunks */
    unsigned long start, stop;
    /* If "x == 0" and "rank == 0", the "count" input parameter is "1",
     * because rank "0" must hold all coefficients. */
    partition_interval_ulong(0, nmax, size, rank, &start, &stop);


    unsigned long *local_order = (unsigned long *)malloc(2 * local_nchunk *
                                                        sizeof(unsigned long));
    if (local_order == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }
    for (size_t j = 0; j < local_nchunk; j++)
        partition_interval_ulong(start, stop, (unsigned long)local_nchunk,
                                 (unsigned long)j, &local_order[2 * j],
                                 &local_order[2 * j + 1]);


    return local_order;
}

