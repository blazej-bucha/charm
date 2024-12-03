/* This header file is not a part of API. */


#ifndef __CHECK_MPI_SHC_H__
#define __CHECK_MPI_SHC_H__


#include <config.h>
#include "../src/prec.h"
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif


#ifdef __cplusplus
extern "C"
{
#endif


extern long int check_mpi_shc(CHARM(shc) *,
                              const char *,
                              CHARM(err) *err);


#ifdef __cplusplus
}
#endif


#endif
