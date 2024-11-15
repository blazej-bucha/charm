/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"
#include "../crd/crd_point_dh2_chunk.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "mpi_crd_point_quad.h"
#include "mpi_err_gather.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(mpi_crd_point_dh2)(unsigned long nmax,
                                       REAL r,
                                       size_t local_nlat,
                                       size_t local_0_start,
                                       MPI_Comm comm,
                                       CHARM(err) *err)
{
    CHARM(point) *pnt = CHARM(mpi_crd_point_quad)(nmax,
                                                  r,
                                                  local_nlat,
                                                  local_0_start,
                                                  comm,
                                                  CHARM(crd_point_dh2_shape),
                                                  CHARM(crd_point_dh2_chunk),
                                                  err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM(mpi_err_isempty)(err))
    {
        CHARM(crd_point_free)(pnt);
        pnt = NULL;
    }


    CHARM(mpi_err_gather)(err);


    return pnt;
}

