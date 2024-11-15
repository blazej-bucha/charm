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
#include "mpi_shc_local2distributed.h"
/* ------------------------------------------------------------------------- */







void CHARM(mpi_shc_local2distributed)(CHARM(shc) *shcs,
                                      MPI_Comm comm)
{
    shcs->distributed = 1;
    shcs->comm        = comm;


    return;
}
