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
/* ------------------------------------------------------------------------- */






/* Returns "1" if "comm1" and "comm2" are identical or congruent MPI
 * communicators and "0" otherwise */
_Bool CHARM(mpi_two_comm_match)(MPI_Comm comm1,
                                MPI_Comm comm2)
{
    int match;
    MPI_Comm_compare(comm1, comm2, &match);


    if ((match == MPI_IDENT) || (match == MPI_CONGRUENT))
        return 1;
    else
        return 0;
}
