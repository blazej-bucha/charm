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
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "mpi_crd_point_check_struct.h"
#include "mpi_crd_point_local2distributed.h"
#include "mpi_err_gather.h"
#include "mpi_err_isdistributed.h"
#include "mpi_crd_point_init_base.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(mpi_crd_point_init_base)(int type,
                                             size_t local_nlat,
                                             size_t local_nlon,
                                             size_t local_0_start,
                                             REAL *lat,
                                             REAL *lon,
                                             REAL *r,
                                             _Bool check_symm,
                                             MPI_Comm comm,
                                             CHARM(err) *err)
{
    /* ----------------------------------------------------------------- */
    CHARM(point) *pnt = NULL;


    if (!CHARM(mpi_err_isdistributed)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }
    /* ----------------------------------------------------------------- */


    /* Prepare the output structure */
    /* ----------------------------------------------------------------- */
    pnt = CHARM(crd_point_init)(type, local_nlat, local_nlon, lat, lon, r);
    if (pnt == NULL)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Failed to create the local portion of the "
                       "\"charm" CHARM_SUFFIX "_point\" structure.");


BARRIER:
    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;


    CHARM(mpi_crd_point_local2distributed)(pnt, local_0_start, comm);
    /* ----------------------------------------------------------------- */


    /* Check the distribution of the latitudinal chunks */
    /* ----------------------------------------------------------------- */
    CHARM(mpi_crd_point_check_struct)(pnt, check_symm, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* ----------------------------------------------------------------- */


EXIT:
    CHARM(mpi_err_gather)(err);


    return pnt;


FAILURE:
    CHARM(crd_point_free)(pnt);
    pnt = NULL;


    goto EXIT;
}
