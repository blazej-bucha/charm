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
#include "testing_msg.h"
#include "test_suite_start.h"
#include "test_suite_end.h"
#include "module_mpi.h"
/* ------------------------------------------------------------------------- */






int main(int argc, char *argv[])
{
#if CHARM_OPENMP
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided != MPI_THREAD_FUNNELED)
    {
        fprintf(stderr, "MPI did not provide MPI_THREAD_FUNNELED\n");
        exit(CHARM_FAILURE);
    }
#else
    MPI_Init(&argc, &argv);
#endif


    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    if (rank == 0)
    {
        test_suite_start();
        printf("Using %d MPI %s\n",
               size, size == 1 ? "process" : "processes");
        printf("\n");


        TESTING_MSG("mpi");
    }


    long int e = module_mpi();


    if (rank == 0)
        test_suite_end(e);


    MPI_Finalize();


    return CHARM_SUCCESS;
}

