/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"
#include "../err/err_propagate.h"
#include "mpi_crd_point_alloc.h"
#include "mpi_err_gather.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(mpi_crd_point_malloc)(int type,
                                          size_t local_nlat,
                                          size_t local_nlon,
                                          size_t local_0_start,
                                          MPI_Comm comm,
                                          CHARM(err) *err)
{
    CHARM(point) *ret = CHARM(mpi_crd_point_alloc)(type, local_nlat,
                                                   local_nlon, local_0_start,
                                                   comm, malloc, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM(mpi_err_isempty)(err))
    {
        CHARM(crd_point_free)(ret);
        ret = NULL;
    }


    CHARM(mpi_err_gather)(err);


    return ret;
}

