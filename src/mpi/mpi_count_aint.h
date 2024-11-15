/* This header file is not a part of API.
 *
 * */


#ifndef __MPI_COUNT_AINT_H__
#define __MPI_COUNT_AINT_H__


#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif


#undef CHARM_MPI_COUNT
#undef CHARM_MPI_AINT
#undef CHARM_MPI_COUNT_DATATYPE
#undef CHARM_MPI_AINT_DATATYPE
#undef USE__C_MPI_ROUTINES
#if HAVE_MPI_IBCAST_C && HAVE_MPI_IREDUCE_C && HAVE_MPI_ALLGATHER_C && \
    HAVE_MPI_GATHERV_V
#   define USE__C_MPI_ROUTINES      1
#   define CHARM_MPI_COUNT          MPI_Count
#   define CHARM_MPI_AINT           MPI_Aint
#   define CHARM_MPI_COUNT_DATATYPE MPI_COUNT
#   define CHARM_MPI_AINT_DATATYPE  MPI_AINT
#else
#   define USE__C_MPI_ROUTINES      0
#   define CHARM_MPI_COUNT          int
#   define CHARM_MPI_AINT           int
#   define CHARM_MPI_COUNT_DATATYPE MPI_INT
#   define CHARM_MPI_AINT_DATATYPE  int
#endif


#endif

