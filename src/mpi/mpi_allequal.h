/* This header file is not a part of API. */


#ifndef __MPI_ALLEQUAL_H__
#define __MPI_ALLEQUAL_H__


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


extern _Bool CHARM(mpi_allequal__Bool)(_Bool,
                                       MPI_Comm,
                                       CHARM(err) *);


extern _Bool CHARM(mpi_allequal_int)(int,
                                     MPI_Comm,
                                     CHARM(err) *);


extern _Bool CHARM(mpi_allequal_ulong)(unsigned long,
                                       MPI_Comm,
                                       CHARM(err) *);


extern _Bool CHARM(mpi_allequal_size_t)(size_t,
                                        MPI_Comm,
                                        CHARM(err) *);


extern _Bool CHARM(mpi_allequal_real)(REAL,
                                      MPI_Comm,
                                      CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
