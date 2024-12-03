/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "check_func.h"
/* ------------------------------------------------------------------------- */






void check_func(const char *func)
{
#if HAVE_MPI
    int rank = 0;
    int initialized;
    MPI_Initialized(&initialized);
    if (initialized)
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (rank == 0)
#endif
    {
        printf("    %s... ", func);
        fflush(stdout);
    }


    return;
}

