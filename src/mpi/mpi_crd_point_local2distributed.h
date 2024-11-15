/* This header file is not a part of API. */


#ifndef __MPI_CRD_POINT_LOCAL2DISTRIBUTED_H__
#define __MPI_CRD_POINT_LOCAL2DISTRIBUTED_H__


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


extern void CHARM(mpi_crd_point_local2distributed)(CHARM(point) *pnt,
                                                   size_t local_0_start,
                                                   MPI_Comm comm);


#ifdef __cplusplus
}
#endif


#endif
