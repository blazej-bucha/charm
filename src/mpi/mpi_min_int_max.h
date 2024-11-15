/* This header file is not a part of API. */


#ifndef __MPI_MIN_INT_MAX_H__
#define __MPI_MIN_INT_MAX_H__


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


extern int CHARM(mpi_min_int_max)(MPI_Comm);


#ifdef __cplusplus
}
#endif


#endif
