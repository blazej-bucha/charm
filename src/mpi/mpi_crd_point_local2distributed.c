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
#include "../crd/crd_point_npoint.h"
#include "../crd/crd_point_isSctr.h"
#include "mpi_size_t.h"
#include "mpi_crd_point_local2distributed.h"
/* ------------------------------------------------------------------------- */







void CHARM(mpi_crd_point_local2distributed)(CHARM(point) *pnt,
                                            size_t local_0_start,
                                            MPI_Comm comm)
{
    /* "pnt" was created locally, which means "pnt->npoint" now stores the
     * total number of points available to an MPI process.  So copy this value
     * to "pnt->local_npoint". */
    pnt->local_npoint  = pnt->npoint;


    /* Now compute the global number of latitudes, longitudes and points */
    size_t local_nlat_sum;
    MPI_Allreduce(&(pnt->local_nlat), &local_nlat_sum, 1, CHARM_MPI_SIZE_T,
                  MPI_SUM, comm);
    pnt->nlat = local_nlat_sum;
    if (CHARM(crd_point_isSctr)(pnt->type))
    {
        size_t local_nlon_sum;
        /* With scattered points, we need to update also "pnt->nlon", because
         * each process may have different number of points
         * ("pnt->local_nlat"). */
        MPI_Allreduce(&(pnt->local_nlon), &local_nlon_sum, 1, CHARM_MPI_SIZE_T,
                      MPI_SUM, comm);
        pnt->nlon = local_nlon_sum;
    }
    pnt->npoint = CHARM(crd_point_npoint)(pnt->type, pnt->nlat, pnt->nlon);


    pnt->local_0_start = local_0_start;
    pnt->distributed   = 1;
    pnt->comm          = comm;


    return;
}
