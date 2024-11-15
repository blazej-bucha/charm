/* This header file is not a part of API. */


#ifndef __MPI_SHC_DISTRIBUTE_H__
#define __MPI_SHC_DISTRIBUTE_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern unsigned long * mpi_shc_distribute(unsigned long,
                                          int,
                                          int,
                                          size_t);


#ifdef __cplusplus
}
#endif


#endif
