/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"
#include "mpi_min_int_max.h"
/* ------------------------------------------------------------------------- */






/* Returns the minimum value of "INT_MAX" across all processes in "comm" */
int CHARM(mpi_min_int_max)(MPI_Comm comm)
{
    int min_INT_MAX;


    int local_INT_MAX = INT_MAX;
    MPI_Allreduce(&local_INT_MAX, &min_INT_MAX, 1, MPI_INT, MPI_MIN, comm);


    return min_INT_MAX;
}
