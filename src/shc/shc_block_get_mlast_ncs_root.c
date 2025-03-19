/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include "../prec.h"
#if HAVE_MPI
#   include <mpi.h>
#   include "../mpi/mpi_real.h"
#   if !HAVE_MPI_BCAST_C
#       include "../mpi/mpi_min_int_max.h"
#   endif
#endif
#include "../err/err_set.h"
#include "../err/err_isempty_all_mpi_processes.h"
#include "shc_block_struct.h"
#include "shc_block_nan.h"
#include "shc_block_get_mlast_ncs_root.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_block_get_mlast_ncs_root)(const CHARM(shc) *shcs,
                                         CHARM(shc_block) *shcs_block,
                                         unsigned long m_get,
                                         unsigned long *mlast,
                                         size_t *ncs,
                                         int *root,
                                         CHARM(err) *err)
{
    char err_msg[CHARM_ERR_MAX_MSG];


    if (m_get > shcs->nmax)
    {
        snprintf(err_msg, CHARM_ERR_MAX_MSG,
                         "\"m_get = %lu\" cannot be larger than "
                         "\"shcs->nmax = %lu\".", m_get, shcs->nmax);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
    }


    if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
        goto FAILURE;


#if HAVE_MPI
    if (shcs->distributed)
    {
        /* ----------------------------------------------------------------- */
        int size;
        MPI_Comm_size(shcs_block->comm, &size);
        /* ----------------------------------------------------------------- */


        /* Identify the MPI process and a its chunk that holds coefficients of
         * order "m_get". */
        /* ----------------------------------------------------------------- */
        _Bool have_m = 0;
        unsigned long chunk_mfirst, chunk_mlast;  /* Minimum and maximum chunk
                                                   * orders */
        for (size_t j = 0; j < shcs->local_nchunk; j++)
        {
            chunk_mfirst = shcs->local_order[2 * j];
            chunk_mlast  = shcs->local_order[2 * j + 1];
            if ((m_get >= chunk_mfirst) && (m_get <= chunk_mlast))
            {
                have_m = 1;
                break;
            }
        }
        /* ----------------------------------------------------------------- */


        /* Gather "have_m" and find out for which process "have_m == 1" */
        /* ----------------------------------------------------------------- */
        MPI_Allgather(&have_m, 1, MPI_C_BOOL, shcs_block->have_m_all, 1,
                      MPI_C_BOOL, shcs_block->comm);


        /* Check that "have_m == 1" on one process only */
        int sum = 0;
        for (int i = 0; i < size; i++)
            sum += (int)shcs_block->have_m_all[i];
        if (sum != 1)
        {
            snprintf(err_msg, CHARM_ERR_MAX_MSG,
                             "\"%d\" processes hold spherical harmonic "
                             "coefficients of order \"%lu\", but this is "
                             "allowed for one process only.", sum, m_get);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
        }


        if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
            goto FAILURE;


        /* Get the rank of the process, on which "have_m == 1" */
        int root_tmp;
        for (root_tmp = 0; root_tmp < size; root_tmp++)
            if (shcs_block->have_m_all[root_tmp])
                break;


        /* "*root_tmp" now holds the rank of the process, for which "have_m ==
         * 1". */
        *root = root_tmp;
        /* ----------------------------------------------------------------- */


        /* ----------------------------------------------------------------- */
#   if !HAVE_MPI_BCAST_C
        /* We can't use "MPI_Bcast_c", so we have to get the minimum value of
         * "INT_MAX" across all processes */
        const int int_max = CHARM(mpi_min_int_max)(shcs_block->comm);


        /* Sanity check */
        if (m_get > (size_t)int_max)
        {
            snprintf(err_msg, CHARM_ERR_MAX_MSG,
                             "\"m_get = %lu\" cannot be larger than "
                             "\"INT_MAX = %d\".", m_get, int_max);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
        }


        if (!CHARM_ERR_ISEMPTY_ALL_MPI_PROCESSES(err))
            goto FAILURE;
#   endif
        /* ----------------------------------------------------------------- */


        /* Get the amount of coefficients that will be sent/received by MPI
         * processes */
        /* ----------------------------------------------------------------- */
        /* Make sure all processes know the maximum order of the coefficients
         * that will be sent/received.  The minimum order is given by the input
         * "m_get" value. */
        MPI_Bcast(&chunk_mlast, 1, MPI_UNSIGNED_LONG, *root, shcs_block->comm);


        /* Get the amount of coefficients to be transmitted */
        unsigned long m;
        size_t ncoeffs     = 0;
        size_t ncoeffs_tmp = 0;
        size_t ncoeffs_max = shcs_block->ncs_max;
        /* The loop runs up to "chunk_mlast", which is the maximum order of the
         * chunk on the root, from which we are getting the coefficients */
        for (m = m_get; m <= chunk_mlast; m++)
        {
            ncoeffs_tmp = ncoeffs + (shcs->nmax - m) + 1;
            if ((ncoeffs_tmp <= ncoeffs_max)
#   if !HAVE_MPI_BCAST_C
                /* If we do not have "MPI_Bcast_c", we have to additionally
                 * make sure that we do not exceed the range of "INT_MAX" */
                && (ncoeffs_tmp <= (size_t)int_max)
#   endif
               )
                ncoeffs = ncoeffs_tmp;
            else
                break;
        }
        m--;


        /* "m" now represents the maximum order of coefficients that will be
         * transmitted.  It might be smaller (but never larger) than
         * "chunk_mlast" because of "CHARM(glob_shc_block_nmax_multiplier)" and
         * "int_max".  "ncoeffs" is the number of "C" or "S" coefficients to be
         * transmitted. */
        *mlast = m;
        *ncs   = ncoeffs;
        /* ----------------------------------------------------------------- */
    }
#endif


    if (!shcs->distributed)
    {
        *mlast = shcs->nmax;
        *ncs   = shcs->nc;  /* "shcs->nc == shcs->ns" */
        *root  = 0;
    }


EXIT:
    return;


FAILURE:
    CHARM(shc_block_nan)(shcs_block);


    *mlast = ULONG_MAX;
    *ncs   = SIZE_MAX;
    *root  = -1;


    goto EXIT;
}
