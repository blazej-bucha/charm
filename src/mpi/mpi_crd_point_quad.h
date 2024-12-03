/* This header file is not a part of API. */


#ifndef __MPI_CRD_POINT_QUAD_H__
#define __MPI_CRD_POINT_QUAD_H__


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


extern CHARM(point) *CHARM(mpi_crd_point_quad)(unsigned long,
                                               REAL,
                                               size_t,
                                               size_t,
                                               MPI_Comm,
                                               void (*)(unsigned long,
                                                        size_t *,
                                                        size_t *),
                                               CHARM(point) *(*)(unsigned long,
                                                                 REAL,
                                                                 size_t,
                                                                 size_t,
                                                                 CHARM(err) *),
                                               CHARM(err) *err);


#ifdef __cplusplus
}
#endif


#endif
