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
#include "mpi_two_comm_match.h"
#include "mpi_three_comm_match.h"
/* ------------------------------------------------------------------------- */






/* Returns "1" if "comm1", "comm2" and "comm3" are identical or congruent MPI
 * communicators and "0" otherwise */
_Bool CHARM(mpi_three_comm_match)(MPI_Comm comm1,
                                  MPI_Comm comm2,
                                  MPI_Comm comm3)
{
    return CHARM(mpi_two_comm_match)(comm1, comm2) &&
           CHARM(mpi_two_comm_match)(comm1, comm3);
}
