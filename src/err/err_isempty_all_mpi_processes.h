/* This header file is not a part of API. */


#ifndef __ERR_ISEMPTY_ALL_MPI_PROCESSES_H__
#define __ERR_ISEMPTY_ALL_MPI_PROCESSES_H__


#include <config.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"


/* "CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES" calls either "err_isempty" or
 * "mpi_err_isempty", depending on whether or not CHarm is being compilied with
 * MPI support.  Importantly, "CHARM(mpi_err_isempty)" is a collective call. */
#undef CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES
#if HAVE_MPI
#   define CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(x) CHARM(mpi_err_isempty)((x))
#else
#   define CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(x) CHARM(err_isempty)((x))
#endif


#endif
