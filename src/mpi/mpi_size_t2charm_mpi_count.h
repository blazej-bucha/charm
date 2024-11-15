/* This header file is not a part of API. */


#ifndef __MPI_SIZE_T2CHARM_MPI_COUNT_H__
#define __MPI_SIZE_T2CHARM_MPI_COUNT_H__


#include <config.h>
#include "../prec.h"
#include "mpi_count_aint.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern CHARM_MPI_COUNT CHARM(mpi_size_t2charm_mpi_count)(size_t,
                                                         CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
