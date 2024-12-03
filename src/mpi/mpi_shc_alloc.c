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
#include "../shc/shc_check_chunk_orders.h"
#include "../shc/shc_local_ncs.h"
#include "../err/err_propagate.h"
#include "../err/err_set.h"
#include "mpi_err_isdistributed.h"
#include "mpi_shc_alloc.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(mpi_shc_alloc)(unsigned long nmax,
                                 REAL mu,
                                 REAL r,
                                 size_t nchunk,
                                 const unsigned long *order,
                                 MPI_Comm comm,
                                 void *(*alloc)(size_t),
                                 CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs = NULL;
    REAL *c          = NULL;
    REAL *s          = NULL;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    if (!CHARM(mpi_err_isdistributed)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }
    /* --------------------------------------------------------------------- */


    /* Total number of coefficients to be stored by an MPI process */
    /* --------------------------------------------------------------------- */
    CHARM(shc_check_chunk_orders)(nmax, nchunk, order, 0, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }


    size_t local_ncs = CHARM(shc_local_ncs)(nmax, nchunk, order);
    /* --------------------------------------------------------------------- */


    /* Allocate the local portion of the coefficients */
    /* --------------------------------------------------------------------- */
    c = (REAL *)alloc(local_ncs * sizeof(REAL));
    if (c == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto BARRIER;
    }


    s = (REAL *)alloc(local_ncs * sizeof(REAL));
    if (s == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto BARRIER;
    }


BARRIER:
    /* Barrier before a collective call */
    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* Create the structure */
    /* --------------------------------------------------------------------- */
    shcs = CHARM(mpi_shc_init)(nmax, mu, r, c, s, nchunk, order, comm, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;


    shcs->owner = 1;
    /* --------------------------------------------------------------------- */


EXIT:
    return shcs;


FAILURE:
    free(c);
    free(s);
    CHARM(shc_free)(shcs);
    shcs = NULL;


    goto EXIT;
}
