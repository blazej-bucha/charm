/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#if HAVE_MPI
#   include <mpi.h>
#endif
#include "../prec.h"
#include "../misc/misc_check_radius.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "shc_check_chunk_orders.h"
#include "shc_local_ncs.h"
#include "shc_init_chunk.h"
/* ------------------------------------------------------------------------- */






CHARM(shc) *CHARM(shc_init_chunk)(unsigned long nmax,
                                  REAL mu,
                                  REAL r,
                                  REAL *c,
                                  REAL *s,
                                  size_t nchunk,
                                  const unsigned long *order,
                                  CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    CHARM(shc) *shcs = (CHARM(shc) *)malloc(sizeof(CHARM(shc)));
    if (shcs == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        return shcs;
    }


    shcs->c = shcs->s = NULL;
#if HAVE_MPI
    shcs->local_order = NULL;
#endif
    /* --------------------------------------------------------------------- */


    /* Save the scalar input values to the "CHARM(shc)" structure */
    /* --------------------------------------------------------------------- */
    shcs->nmax = nmax;
    shcs->mu   = mu;


    CHARM(misc_check_radius)(r, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto FAILURE;
    }
    shcs->r = r;


    /* Size of the first dimension of the "shcs->c" and "shcs->s" pointer
     * arrays (number of harmonic orders) */
    size_t nmax1 = (size_t)(nmax) + 1;
    /* Total number of spherical harmonic coefficients to be stored in both
     * "shcs->c[0]" and "shcs->s[0]".  Note that the "S_{n,0}" coefficients,
     * which do not exist, will be stored in "shcs->s[0]" as zeros. */
    size_t len = ((nmax1 + 1) * nmax1) / 2;
    shcs->nc = shcs->ns = len;


    shcs->owner       = 0;
    shcs->distributed = 0;
    /* --------------------------------------------------------------------- */


    /* Allocate the "shcs->c" and "shcs->s" pointer arrays */
    /* --------------------------------------------------------------------- */
    shcs->c = (REAL **)malloc(nmax1 * sizeof(REAL *));
    if (shcs->c == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }


    shcs->s = (REAL **)malloc(nmax1 * sizeof(REAL *));
    if (shcs->s == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }


    for (size_t m = 0; m <= nmax; m++)
        shcs->c[m] = shcs->s[m] = NULL;
    /* --------------------------------------------------------------------- */


    /* If compiling with the MPI support, return structure with no coefficients
     * if "nchunk" is "0".  If compiling without the MPI support, report
     * a failure if "nchunk == 0". */
    /* --------------------------------------------------------------------- */
    if (nchunk == 0)
    {
#if HAVE_MPI
        shcs->local_nc = shcs->local_ns = 0;
        shcs->local_nchunk = nchunk;
        shcs->local_order = NULL;


        goto EXIT;
#else
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"nchunk\" must be larger than zero.");
        goto FAILURE;
#endif
    }
    /* --------------------------------------------------------------------- */


    /* Set "shcs->c[m]" and "shcs->s[m]" to point to "C_{m,m}" and "S_{m,m}" */
    /* --------------------------------------------------------------------- */
    CHARM(shc_check_chunk_orders)(nmax, nchunk, order, 0, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto FAILURE;
    }


    /* The "c" and "s" pointers of the input arrays are always attached to
     * "shcs->c[order[0]]" and "shcs->s[order[0]]", respectively.  This is
     * important when implementing "shc_free". */
    size_t idx = 0;
    for (size_t j = 0; j < nchunk; j++)
    {
        for (unsigned long m = order[2 * j]; m <= order[2 * j + 1]; m++)
        {
            shcs->c[m] = c + idx;
            shcs->s[m] = s + idx;


            idx += nmax1 - m;
        }
    }
    /* --------------------------------------------------------------------- */


    /* Set the MPI-related members of "shcs" */
    /* --------------------------------------------------------------------- */
#if HAVE_MPI
    shcs->local_nc = shcs->local_ns = CHARM(shc_local_ncs)(nmax, nchunk,
                                                           order);
    shcs->local_nchunk = nchunk;


    unsigned long *local_order = (unsigned long *)
                                    malloc(2 * nchunk * sizeof(unsigned long));
    if (local_order == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto FAILURE;
    }
    for (size_t j = 0; j < 2 * nchunk; j++)
        local_order[j] = order[j];


    shcs->local_order = local_order;
    shcs->comm        = MPI_COMM_NULL;
#endif
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
EXIT:
    return shcs;
    /* --------------------------------------------------------------------- */


FAILURE:
    /* --------------------------------------------------------------------- */
    free(shcs->c);
    free(shcs->s);
#if HAVE_MPI
    free(shcs->local_order);
#endif
    free(shcs);
    shcs = NULL;


    goto EXIT;
    /* --------------------------------------------------------------------- */
}
