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
#include "mpi_two_comm_match.h"
#include "mpi_check_point_shc_err.h"
/* ------------------------------------------------------------------------- */






void CHARM(mpi_check_point_shc_err)(const CHARM(point) *pnt,
                                    const CHARM(shc) *shcs,
                                    CHARM(err) *err)
{
    if (pnt->distributed != shcs->distributed)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"pnt\" and \"shcs\" must either be both distributed "
                       "or both non-distributed.");
        return;
    }


    if (pnt->distributed && shcs->distributed && !err->distributed)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"err\" must be distributed if "
                       "\"pnt\" and \"shcs\" are distributed.");
        return;
    }


    if (!pnt->distributed && !shcs->distributed && err->distributed)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"err\" must not be distributed if "
                       "\"pnt\" and \"shcs\" are not distributed ");
        return;
    }


    if (pnt->distributed)
    {
        if (!CHARM(mpi_two_comm_match)(pnt->comm, err->comm))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The MPI communicators \"pnt->comm\" "
                           "and \"err->comm\" must be "
                           "identical or congruent (\"MPI_IDENT\" or "
                           "\"MPI_CONGRUENT\").");
            return;
        }
    }


    if (shcs->distributed)
    {
        if (!CHARM(mpi_two_comm_match)(shcs->comm, err->comm))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The MPI communicators \"shcs->comm\" "
                           "and \"err->comm\" must be "
                           "identical or congruent (\"MPI_IDENT\" or "
                           "\"MPI_CONGRUENT\").");
            return;
        }
    }


    return;
}
