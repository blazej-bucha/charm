/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../crd/crd_point_isDHGrid.h"
#include "../crd/crd_point_isQuadGrid.h"
#include "../crd/crd_point_isCustGrid.h"
#include "../crd/crd_point_isGrid.h"
#include "../crd/crd_point_isSctr.h"
#include "mpi_size_t.h"
#include "mpi_allequal.h"
#include "mpi_err_gather.h"
#include "mpi_allequal.h"
#include "mpi_crd_point_issymm.h"
#include "mpi_crd_point_check_struct.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
#undef ALLEQUAL
#define ALLEQUAL(type, x, xstr)                                               \
    {                                                                         \
        _Bool allequal = CHARM(CAT(mpi_allequal_, type))(x, pnt->comm, err);  \
        if (!CHARM(err_isempty)(err))                                         \
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);          \
                                                                              \
                                                                              \
        if (!CHARM(mpi_err_isempty)(err))                                     \
            goto EXIT;                                                        \
                                                                              \
                                                                              \
        if (!allequal)                                                        \
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG, \
                           "\"" xstr "\" must be equal for all "              \
                           "processes in \"pnt->comm\".");                    \
                                                                              \
                                                                              \
        if (!CHARM(mpi_err_isempty)(err))                                     \
            goto EXIT;                                                        \
    }
/* ------------------------------------------------------------------------- */






void CHARM(mpi_crd_point_check_struct)(const CHARM(point) *pnt,
                                       _Bool check_symm,
                                       CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    int rank, size;
    MPI_Comm_rank(pnt->comm, &rank);
    MPI_Comm_size(pnt->comm, &size);


    char err_msg[CHARM_ERR_MAX_MSG];


    size_t *local_nlat    = NULL;
    size_t *local_0_start = NULL;
    /* --------------------------------------------------------------------- */


    /* Some values must be equal across all processes */
    /* --------------------------------------------------------------------- */
    ALLEQUAL(int,    pnt->type,        "pnt->type");
    ALLEQUAL(size_t, pnt->nlat,        "pnt->nlat");
    ALLEQUAL(size_t, pnt->nlon,        "pnt->nlon");
    ALLEQUAL(size_t, pnt->npoint,      "pnt->npoint");
    ALLEQUAL(_Bool,  pnt->owner,       "pnt->owner");
    ALLEQUAL(_Bool,  pnt->distributed, "pnt->distributed");
    /* --------------------------------------------------------------------- */


    /* Check "local_nlat" and "local_0_start" */
    /* --------------------------------------------------------------------- */
    local_nlat = (size_t *)malloc(size * sizeof(size_t));
    if (local_nlat == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto BARRIER;
    }


    local_0_start = (size_t *)malloc(size * sizeof(size_t));
    if (local_0_start == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto BARRIER;
    }


    /* Barrier before collective calls */
BARRIER:
    if (!CHARM(mpi_err_isempty)(err))
        goto EXIT;


    MPI_Gather(&(pnt->local_nlat), 1, CHARM_MPI_SIZE_T,
               local_nlat, 1, CHARM_MPI_SIZE_T, 0, pnt->comm);
    MPI_Gather(&(pnt->local_0_start), 1, CHARM_MPI_SIZE_T,
               local_0_start, 1, CHARM_MPI_SIZE_T, 0, pnt->comm);
    _Bool symm = (check_symm) ? CHARM(mpi_crd_point_issymm)(pnt) : 0;


    if (rank == 0)
    {
        /* Check that the sum of "pnt->local_nlat" across all processes is
         * equal to "pnt->nlat".  This could be done with "MPI_Allreduce", but
         * we need "local_nlat" later anyway. */
        /* ----------------------------------------------------------------- */
        {
        size_t sum = 0;
        for (int i = 0; i < size; i++)
            sum += local_nlat[i];
        if (sum != pnt->nlat)
        {
            snprintf(err_msg, CHARM_ERR_MAX_MSG,
                             "The sum of \"pnt->local_nlat\" is \"%zu\" which "
                             "does not match the total number in "
                             "\"pnt->nlat = %zu\".",
                             sum, pnt->nlat);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            goto EXIT;
        }
        }
        /* ----------------------------------------------------------------- */


        /* Check that across all processes with "local_nlat > 0", there is only
         * one process with "pnt->local_0_start == 0".  Processes with
         * "local_nlat == 0" are allowed to have any "local_0_start", even
         * zero; this signalizes that that particular process stores zero
         * points. */
        /* ----------------------------------------------------------------- */
        {
        int sum = 0;
        for (int i = 0; i < size; i++)
            if (local_nlat[i] > 0)
                sum += (int)(local_0_start[i] == 0);


        if (sum > 1)
        {
            snprintf(err_msg, CHARM_ERR_MAX_MSG,
                             "\"pnt->local_0_start\" is \"0\" at \"%d\" "
                             "processes, but this must be true at one "
                             "process.", sum);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            goto EXIT;
        }
        }
        /* ----------------------------------------------------------------- */


        /* Check that the local length of latitudinal chunks matches the global
         * indexing.  For instance,
         *
         *          size_t local_0_start[4] = {0, 3,  8, 12};
         *          size_t local_nlat[4]    = {3, 5,  4,  5};
         *
         * is valid, but
         *
         *          size_t local_0_start[4] = {0, 3,  8, 12};
         *          size_t local_nlat[4]    = {3, 5,  5,  5};
         *
         * is not valid, because the third chunk starting at the global index
         * "8" has size "5", implying the next chunk should start at "13", but
         * in reality it starts at "12". */
        /* ----------------------------------------------------------------- */
        size_t start = 0;
        size_t add = 0;
        _Bool match;


        /* Quadrature and symmetric grids treat "local_nlat" and
         * "local_0_start" differently.  Specifically, "local_nlat" has to be
         * divided by "2" in the algorithm that follows, because we check in
         * this case the northern hemisphere only, hence the number of points
         * in "local_nlat" is halved.  This treats correctly also symmetric
         * grids having the equator. */
        size_t div_factor;
        if (CHARM(crd_point_isQuadGrid)(pnt->type))
            div_factor = 2;
        else if (CHARM(crd_point_isCustGrid)(pnt->type) && symm)
            div_factor = 2;
        else
            div_factor = 1;


        for (int i = 0; i < size; i++)
        {
            if (local_nlat[i] == 0)
                /* This indicates empty structure */
                continue;


            match = 0;
            for (int j = 0; j < size; j++)
            {
                if (local_0_start[j] == start)
                {
                    add = local_nlat[j] / div_factor;
                    if (CHARM(crd_point_isDHGrid)(pnt->type) &&
                        (local_0_start[j] == 0))
                        /* For the first chunk of Driscoll--Healy grids, we
                         * must increase "add" by "1" */
                        add += 1;
                    start += add;
                    match = 1;
                    break;
                }
            }

            if (!match)
            {
                snprintf(err_msg, CHARM_ERR_MAX_MSG,
                                 "Wrong distribution of latitudinal chunks.  "
                                 "Expected to find "
                                 "\"pnt->local_0_start = %zu\" on "
                                 "one MPI process.", start);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG, err_msg);
                goto EXIT;
            }
        }
        /* ----------------------------------------------------------------- */
    }
    /* --------------------------------------------------------------------- */


    /* Check "local_nlon" */
    /* --------------------------------------------------------------------- */
    if (CHARM(crd_point_isGrid)(pnt->type))
    {
        _Bool ret = CHARM(mpi_allequal_size_t)(pnt->local_nlon, pnt->comm,
                                               err);
        if (!CHARM(err_isempty)(err))
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


        if (!CHARM(mpi_err_isempty)(err))
                goto EXIT;


        if (!ret)
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "For all point grids, \"pnt->local_nlon\" must be "
                           "equal for all processes in \"pnt->comm\".");


        if (!CHARM(mpi_err_isempty)(err))
            goto EXIT;
    }
    else if (CHARM(crd_point_isSctr)(pnt->type))
    {
        size_t local_nlon_sum;
        MPI_Allreduce(&(pnt->local_nlon), &local_nlon_sum, 1, CHARM_MPI_SIZE_T,
                      MPI_SUM, pnt->comm);
        if (local_nlon_sum != pnt->nlon)
        {
            snprintf(err_msg, CHARM_ERR_MAX_MSG,
                             "The sum of \"pnt->local_nlon = %zu\" across "
                             "all processes in \"pnt->comm\" does not match "
                             "the total number of latitudes "
                             "\"pnt->nlon = %zu\".",
                             local_nlon_sum, pnt->nlon);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                           CHARM_EFUNCARG, err_msg);


            if (!CHARM(mpi_err_isempty)(err))
                goto EXIT;
        }
    }
    /* --------------------------------------------------------------------- */


EXIT:
    free(local_nlat);
    free(local_0_start);


    return;
}
