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
#include "mpi_size_t.h"
#include "mpi_allequal.h"
#include "mpi_shc_check_struct.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
#undef ALLEQUAL
#define ALLEQUAL(type, x, xstr)                                               \
    {                                                                         \
        _Bool allequal = CHARM(CAT(mpi_allequal_, type))(x, shcs->comm, err); \
        if (!CHARM(err_isempty)(err))                                         \
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);          \
                                                                              \
                                                                              \
        if (!CHARM(mpi_err_isempty)(err))                                     \
            goto EXIT;                                                        \
                                                                              \
                                                                              \
        if (!allequal)                                                        \
        {                                                                     \
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG, \
                           "\"" xstr "\" must be equal for all "              \
                           "processes in \"shcs->comm\".");                   \
            goto EXIT;                                                        \
        }                                                                     \
    }
/* ------------------------------------------------------------------------- */






void CHARM(mpi_shc_check_struct)(const CHARM(shc) *shcs,
                                 CHARM(err) *err)
{
    int rank, size;
    MPI_Comm_rank(shcs->comm, &rank);
    MPI_Comm_size(shcs->comm, &size);


    char err_msg[CHARM_ERR_MAX_MSG];


    /* Some values must be equal across all processes */
    /* --------------------------------------------------------------------- */
    ALLEQUAL(ulong,  shcs->nmax,        "shcs->nmax");
    ALLEQUAL(real,   shcs->mu,          "shcs->mu");
    ALLEQUAL(real,   shcs->r,           "shcs->r");
    ALLEQUAL(size_t, shcs->nc,          "shcs->nc");
    ALLEQUAL(size_t, shcs->ns,          "shcs->ns");
    ALLEQUAL(_Bool,  shcs->owner,       "shcs->owner");
    ALLEQUAL(_Bool,  shcs->distributed, "shcs->distributed");
    /* --------------------------------------------------------------------- */



    /* Some other values must yield specific sums */
    /* --------------------------------------------------------------------- */
    size_t local_nc, local_ns, local_nchunk;


    MPI_Allreduce(&(shcs->local_nc), &local_nc, 1, CHARM_MPI_SIZE_T, MPI_SUM,
                  shcs->comm);
    MPI_Allreduce(&(shcs->local_ns), &local_ns, 1, CHARM_MPI_SIZE_T, MPI_SUM,
                  shcs->comm);
    MPI_Allreduce(&(shcs->local_nchunk), &local_nchunk, 1, CHARM_MPI_SIZE_T,
                  MPI_SUM, shcs->comm);


    if (local_nc != shcs->nc)
    {
        sprintf(err_msg, "The sum of \"shcs->local_nc\" across MPI "
                         "processes in \"shcs->comm\" is \"%zu\" which "
                         "does not match the total number in "
                         "\"shcs->nc = %zu\".",
                         local_nc, shcs->nc);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto EXIT;
    }


    if (local_ns != shcs->ns)
    {
        sprintf(err_msg, "The sum of \"shcs->local_ns\" across MPI "
                         "processes in \"shcs->comm\" is \"%zu\" which "
                         "does not match the total number in "
                         "\"shcs->ns = %zu\".",
                         local_ns, shcs->ns);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto EXIT;
    }


    if (local_nchunk == 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The sum of \"shcs->local_nchunk\" across MPI "
                       "processes in \"shcs->comm\" cannot be zero.");
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */


EXIT:
    return;
}
