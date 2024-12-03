/* This header file is not a part of API. */


#ifndef __MPI_THREE_COMM_MATCH_H__
#define __MPI_THREE_COMM_MATCH_H__


#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern _Bool CHARM(mpi_three_comm_match)(MPI_Comm,
                                         MPI_Comm,
                                         MPI_Comm);


#ifdef __cplusplus
}
#endif


#endif
