/* This header file is not a part of API.
 *
 * */


#ifndef __MPI_ALLGATHER_H__
#define __MPI_ALLGATHER_H__


#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "mpi_count_aint.h"


#undef CHARM_MPI_ALLGATHER
/* We required both "HAVE_MPI_BCAST_C" and "HAVE_MPI_REDUCE_C" to simplify
 * things because of "mpi_count_aint.h" */
#if USE__C_MPI_ROUTINES
#    define CHARM_MPI_ALLGATHER  MPI_Allgather_c
#else
#    define CHARM_MPI_ALLGATHER  MPI_Allgather
#endif


#endif

