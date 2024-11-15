/* This header file is not a part of API. */


#ifndef __MPI_CRD_POINT_ALLOC_H__
#define __MPI_CRD_POINT_ALLOC_H__


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


extern CHARM(point) *CHARM(mpi_crd_point_alloc)(int type,
                                                size_t local_nlat,
                                                size_t nlon,
                                                size_t local_0_start,
                                                MPI_Comm comm,
                                                void *(*alloc)(size_t),
                                                CHARM(err) *err);


#ifdef __cplusplus
}
#endif


#endif
