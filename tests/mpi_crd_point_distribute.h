/* This header file is not a part of API. */


#ifndef __MPI_CRD_POINT_DISTRIBUTE_H__
#define __MPI_CRD_POINT_DISTRIBUTE_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void mpi_crd_point_distribute(size_t,
                                     int,
                                     int,
                                     int,
                                     size_t *,
                                     size_t *);


#ifdef __cplusplus
}
#endif


#endif
