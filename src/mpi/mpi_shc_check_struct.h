/* This header file is not a part of API. */


#ifndef __MPI_SHC_CHECK_STRUCT_H__
#define __MPI_SHC_CHECK_STRUCT_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(mpi_shc_check_struct)(const CHARM(shc) *,
                                        CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
