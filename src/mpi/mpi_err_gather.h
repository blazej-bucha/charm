/* This header file is not a part of API. */


#ifndef __MPI_ERR_GATHER_H__
#define __MPI_ERR_GATHER_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(mpi_err_gather)(CHARM(err) *err);


#undef ERR_GATHER_MSG
#define ERR_GATHER_MSG "An error was reported by one or more processes " \
                       "other than this one."


#ifdef __cplusplus
}
#endif


#endif
