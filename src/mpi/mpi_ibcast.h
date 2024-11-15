/* This header file is not a part of API.
 *
 * */


#ifndef __MPI_IBCAST_H__
#define __MPI_IBCAST_H__


#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "mpi_count_aint.h"


#undef CHARM_MPI_IBCAST
#if USE__C_MPI_ROUTINES
#    define CHARM_MPI_IBCAST  MPI_Ibcast_c
#else
#    define CHARM_MPI_IBCAST  MPI_Ibcast
#endif


#endif

