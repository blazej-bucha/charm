/* This header file is not a part of API.
 *
 * */


#ifndef __MPI_REAL_H__
#define __MPI_REAL_H__


#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"


#undef CHARM_MPI_REAL
#if CHARM_FLOAT
#   define CHARM_MPI_REAL MPI_FLOAT
#elif CHARM_QUAD
#   error "The MPI standard does not define __float128 data type"
#else
#   define CHARM_MPI_REAL MPI_DOUBLE
#endif


#endif

