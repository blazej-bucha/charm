/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#if HAVE_MPI
#   include <mpi.h>
#   include "../mpi/mpi_allequal.h"
#endif
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_is_null_ptr.h"
#include "../err/err_isempty_all_mpi_processes.h"
#include "../glob/glob_get_shc_block_nmax_multiplier.h"
#include "shc_block_free.h"
#include "shc_block_get_mlast_ncs_root.h"
#include "shc_block_init.h"
/* ------------------------------------------------------------------------- */






CHARM(shc_block) *CHARM(shc_block_init)(const CHARM(shc) *shcs)
{
    /* --------------------------------------------------------------------- */
    CHARM(shc_block) *shcs_block = NULL;
    CHARM(err) *err              = NULL;
    /* --------------------------------------------------------------------- */


    /* Initialize an error structure */
    /* --------------------------------------------------------------------- */
    err = CHARM(err_init)();
    if (CHARM(err_is_null_ptr)(err
#if HAVE_MPI
                               , shcs->distributed,
                               shcs->comm
#endif
                              ))
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    shcs_block = (CHARM(shc_block) *)malloc(sizeof(CHARM(shc_block)));
    if (shcs_block == NULL)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);


    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto FAILURE;


    shcs_block->nmax = shcs->nmax;
    shcs_block->c = shcs_block->s = NULL;
    /* --------------------------------------------------------------------- */


#if HAVE_MPI
    shcs_block->have_m_all  = NULL;
    shcs_block->distributed = shcs->distributed;
    shcs_block->comm        = MPI_COMM_NULL;


    if (shcs->distributed)
    {
        /* "n" is the maximum number of coefficients that can be stored in
         * "shcs_block->c" and "shcs_block->s" */
        size_t n = CHARM(glob_get_shc_block_nmax_multiplier)() *
                   (shcs->nmax + 1);


        shcs_block->c = (REAL *)calloc(n, sizeof(REAL));
        if (shcs_block->c == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);
            goto BARRIER;
        }


        shcs_block->s = (REAL *)calloc(n, sizeof(REAL));
        if (shcs_block->s == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);
            goto BARRIER;
        }


BARRIER:
        if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
            goto FAILURE;


        shcs_block->ncs_max = n;
        shcs_block->owner   = 1;


        int mpi_size;
        MPI_Comm_size(shcs->comm, &mpi_size);
        shcs_block->have_m_all = (_Bool *)malloc(mpi_size * sizeof(_Bool));
        if (shcs_block->have_m_all == NULL)
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                           CHARM_ERR_MALLOC_FAILURE);


        if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
            goto FAILURE;


        shcs_block->comm = shcs->comm;
    }
#endif


    /* Non-distributed "charm_shc" structure regardless of the MPI support. */
    if (!shcs->distributed)
    {
        shcs_block->c       = shcs->c[0];
        shcs_block->s       = shcs->s[0];
        shcs_block->ncs_max = shcs->nc;  /* "shcs->nc == shcs->ns" */
        shcs_block->owner   = 0;
    }


    /* "mfirst", "mlast", "nc", "ns" */
    /* --------------------------------------------------------------------- */
    shcs_block->mfirst = 0;
    unsigned long mlast;
    size_t ncs;
    int root;
    CHARM(shc_block_get_mlast_ncs_root)(shcs, shcs_block, 0, &mlast, &ncs,
                                        &root, err);
    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto FAILURE;


    shcs_block->mlast = mlast;
    shcs_block->nc    = ncs;
    shcs_block->ns    = ncs;
    shcs_block->root  = root;
    /* --------------------------------------------------------------------- */


EXIT:
    CHARM(err_free)(err);
    return shcs_block;


FAILURE:
    CHARM(shc_block_free)(shcs_block);
    shcs_block = NULL;


    goto EXIT;
}
