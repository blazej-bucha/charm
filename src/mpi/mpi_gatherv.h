/* This header file is not a part of API.
 *
 * */


#ifndef __MPI_GATHERV_H__
#define __MPI_GATHERV_H__


#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "mpi_count_aint.h"


#undef CHARM_MPI_GATHERV
/* We required both "HAVE_MPI_BCAST_C" and "HAVE_MPI_REDUCE_C" to simplify
 * things because of "mpi_count_aint.h" */
#if USE__C_MPI_ROUTINES
#    define CHARM_MPI_GATHERV  MPI_Gatherv_c
#else
#    define CHARM_MPI_GATHERV  MPI_Gatherv
#endif


#endif

