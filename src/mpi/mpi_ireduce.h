/* This header file is not a part of API.
 *
 * */


#ifndef __MPI_IREDUCE_H__
#define __MPI_IREDUCE_H__


#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "mpi_count_aint.h"


#undef CHARM_MPI_IREDUCE
#if USE__C_MPI_ROUTINES
#    define CHARM_MPI_IREDUCE MPI_Ireduce_c
#else
#    define CHARM_MPI_IREDUCE MPI_Ireduce
#endif


#endif

