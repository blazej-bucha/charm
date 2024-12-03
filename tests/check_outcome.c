/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include <stdio.h>
#include "check_outcome.h"
/* ------------------------------------------------------------------------- */






void check_outcome(long int err)
{
#if HAVE_MPI
    int rank = 0;
    int initialized;
    MPI_Initialized(&initialized);
    if (initialized)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        long int err_sum;
        MPI_Allreduce(&err, &err_sum, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
        err = err_sum;
    }


    if (rank == 0)
#endif
    {
        if (err)
            printf("\n            %ld possibly wrong result%s\n\n",
                   err, (err == 1) ? "" : "s");
        else
            printf("ok\n");
    }
}

