/* This header file is not a part of API. */


#ifndef __CHECK_MPI_SHC_ALLOC_H__
#define __CHECK_MPI_SHC_ALLOC_H__


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


extern long int check_mpi_shc_alloc(CHARM(shc) *(*)(unsigned long,
                                                    REAL,
                                                    REAL,
                                                    size_t,
                                                    const unsigned long *,
                                                    MPI_Comm,
                                                    CHARM(err) *));


#ifdef __cplusplus
}
#endif


#endif
