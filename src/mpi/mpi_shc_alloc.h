/* This header file is not a part of API. */


#ifndef __MPI_SHC_ALLOC_H__
#define __MPI_SHC_ALLOC_H__


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


extern CHARM(shc) *CHARM(mpi_shc_alloc)(unsigned long,
                                        REAL,
                                        REAL,
                                        size_t,
                                        const unsigned long *,
                                        MPI_Comm,
                                        void *(*)(size_t),
                                        CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
