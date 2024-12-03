/* This header file is not a part of API. */


#ifndef __MPI_SHC_LOCAL2DISTRIBUTED_H__
#define __MPI_SHC_LOCAL2DISTRIBUTED_H__


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


extern void CHARM(mpi_shc_local2distributed)(CHARM(shc) *shcs,
                                             MPI_Comm comm);


#ifdef __cplusplus
}
#endif


#endif
