/* This header file is not a part of API. */


#ifndef __ERR_OMP_MPI_H__
#define __ERR_OMP_MPI_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern int CHARM(err_omp_mpi)(int *,
                              int *,
                              const char *,
                              int err_code,
                              CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
