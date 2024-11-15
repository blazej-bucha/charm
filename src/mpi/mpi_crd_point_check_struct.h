/* This header file is not a part of API. */


#ifndef __MPI_CRD_POINT_CHECK_STRUCT_H__
#define __MPI_CRD_POINT_CHECK_STRUCT_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(mpi_crd_point_check_struct)(const CHARM(point) *,
                                              _Bool,
                                              CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
