/* This header file is not a part of API. */


#ifndef __MPI_CHECK_POINT_SHC_ERR_H__
#define __MPI_CHECK_POINT_SHC_ERR_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(mpi_check_point_shc_err)(const CHARM(point) *,
                                           const CHARM(shc) *,
                                           CHARM(err) *err);


#ifdef __cplusplus
}
#endif


#endif
