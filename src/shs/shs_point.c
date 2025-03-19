/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shs_point_grd.h"
#include "shs_point_sctr.h"
#include "shs_point_gradn.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_isempty_all_mpi_processes.h"
#include "../crd/crd_point_isSctr.h"
#include "../crd/crd_point_isGrid.h"
#if HAVE_MPI
#   include "../mpi/mpi_check_point_shc_err.h"
#   include "../mpi/mpi_err_gather.h"
#endif
/* ------------------------------------------------------------------------- */






void CHARM(shs_point)(const CHARM(point) *pnt,
                      const CHARM(shc) *shcs,
                      unsigned long nmax,
                      REAL *f,
                      CHARM(err) *err)
{
    char err_msg[CHARM_ERR_MAX_MSG];


    /* Some trivial initial error checks */
    /* --------------------------------------------------------------------- */
    if (nmax > shcs->nmax)
    {
        snprintf(err_msg, CHARM_ERR_MAX_MSG,
                         "Maximum harmonic degree of the synthesis "
                         "\"nmax = %lu\" cannot be larger than maximum "
                         "harmonic degree of spherical harmonic "
                         "coefficients \"shcs->nmax = %lu\".",
                         nmax, shcs->nmax);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto BARRIER;
    }


    if (!CHARM(crd_point_isSctr)(pnt->type) &&
        !CHARM(crd_point_isGrid)(pnt->type))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported \"pnt->type\" for spherical harmonic "
                       "synthesis of point values.");
        goto BARRIER;
    }


#if HAVE_MPI
    CHARM(mpi_check_point_shc_err)(pnt, shcs, err);
    if (!CHARM(mpi_err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }
#endif
    /* --------------------------------------------------------------------- */






    /* Do nothing if the total number of points in "pnt" is zero, which is
     * a valid case */
    /* --------------------------------------------------------------------- */
    if (pnt->npoint == 0)
        goto EXIT;
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
BARRIER:
    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto EXIT;
    /* --------------------------------------------------------------------- */






    /* Now do the synthesis */
    /* --------------------------------------------------------------------- */
    if (CHARM(crd_point_isSctr)(pnt->type))
        CHARM(shs_point_sctr)(pnt, shcs, nmax, GRAD_0, GRAD_0, GRAD_0, &f,
                              err);
    else if (CHARM(crd_point_isGrid)(pnt->type))
        CHARM(shs_point_grd)(pnt, shcs, nmax, GRAD_0, GRAD_0, GRAD_0, &f, err);


    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
EXIT:
#if HAVE_MPI
    CHARM(mpi_err_gather)(err);
#endif


    return;
    /* --------------------------------------------------------------------- */
}
