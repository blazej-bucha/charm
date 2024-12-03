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
#include "../crd/crd_point_alloc.h"
#include "mpi_err_isdistributed.h"
#include "mpi_crd_point_init_base.h"
#include "mpi_crd_point_alloc.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(mpi_crd_point_alloc)(int type,
                                         size_t local_nlat,
                                         size_t local_nlon,
                                         size_t local_0_start,
                                         MPI_Comm comm,
                                         void *(*alloc)(size_t),
                                         CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    CHARM(point) *tmp = NULL;
    CHARM(point) *ret = NULL;


    if (!CHARM(mpi_err_isdistributed)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    tmp = CHARM(crd_point_alloc)(type, local_nlat, local_nlon, alloc);
    if (tmp == NULL)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Failed to create the local portion of the "
                       "\"charm" CHARM_SUFFIX "_point\" structure.");


BARRIER:
    /* Barrier before a collective call */
    if (!CHARM(mpi_err_isempty)(err))
        goto EXIT;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    ret = CHARM(mpi_crd_point_init_base)(tmp->type,
                                         tmp->nlat,
                                         tmp->nlon,
                                         local_0_start,
                                         tmp->lat,
                                         tmp->lon,
                                         tmp->r,
                                         0,
                                         comm,
                                         err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM(mpi_err_isempty)(err))
    {
        CHARM(crd_point_free)(ret);
        ret = NULL;
        goto EXIT;
    }


    ret->w = tmp->w;
    ret->owner = 1;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
EXIT:
    /* Free "tmp", but without the allocated "tmp->lat", "tmp->lon", "tmp->r"
     * and "tmp->w" arrays, which are attached to "ret" */
    if (tmp != NULL)
        tmp->owner = 0;
    CHARM(crd_point_free)(tmp);


    return ret;
    /* --------------------------------------------------------------------- */
}
