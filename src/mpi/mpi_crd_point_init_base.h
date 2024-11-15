/* This header file is not a part of API. */


#ifndef __MPI_CRD_POINT_INIT_BASE_H__
#define __MPI_CRD_POINT_INIT_BASE_H__


#include <config.h>
#include "../prec.h"
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif


#ifdef __cplusplus
extern "C"
{
#endif


extern CHARM(point) *CHARM(mpi_crd_point_init_base)(int,
                                                    size_t,
                                                    size_t,
                                                    size_t,
                                                    REAL *,
                                                    REAL *,
                                                    REAL *,
                                                    _Bool,
                                                    MPI_Comm,
                                                    CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
