/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "shc_check_chunk_orders.h"
/* ------------------------------------------------------------------------- */






/* Checks whether the "order" array properly specifies "nchunk" chunks of
 * spherical harmonic coefficients up to degree "nmax".
 *
 * If "check_gaps == 1", checked is also that there are no gaps between the
 * chunks.  This is useful with the MPI support.  First, one sums "nchunk" and
 * gathers all "order" across all MPI processes.  Then, one uses this routine
 * with "check_gaps = 1" to check for gaps with these news "nchunk_all" and
 * "order_all".
 *
 * */
int CHARM(shc_check_chunk_orders)(unsigned long nmax,
                                  size_t nchunk,
                                  const unsigned long *order,
                                  _Bool check_gaps,
                                  CHARM(err) *err)
{
    char err_msg[CHARM_ERR_MAX_MSG];


    /* No value of "order" can be larger than "nmax" */
    for (size_t j = 0; j < 2 * nchunk; j++)
    {
        if (order[j] > nmax)
        {
            sprintf(err_msg, "Chunk order \"%lu\" is larger than "
                             "\"nmax = %lu\".", order[j], nmax);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            return -1;
        }
    }


    /* The minimum chunk order cannot be larger than the maximum chunk order */
    unsigned long min, max;
    for (size_t j = 0; j < nchunk; j++)
    {
        min = order[2 * j];
        max = order[2 * j + 1];


        if (min > max)
        {
            sprintf(err_msg, "The minimum chunk order \"%lu\" is larger "
                             "than its corresponding maximum chunk order "
                             "\"%lu\".", min, max);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            return -2;
        }
    }


    /* No two chunks can share the same value of minimum or maximum order */
    unsigned long min_tmp, max_tmp;
    for (size_t j = 0; j < nchunk; j++)
    {
        min = order[2 * j];
        max = order[2 * j + 1];


        for (size_t jj = j + 1; jj < nchunk; jj++)
        {
            min_tmp = order[2 * jj];
            max_tmp = order[2 * jj + 1];
            if ((min == min_tmp) || (max == max_tmp))
            {
                sprintf(err_msg, "Chunks \"%lu, %lu\" and \"%lu, %lu\" share "
                                 "the same value of the minimum or maximum "
                                 "order.", min, max, min_tmp, max_tmp);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG, err_msg);
                return -3;
            }
        }
    }


    /* The chunks must not overlap */
    for (size_t j = 0; j < nchunk; j++)
    {
        min = order[2 * j];
        max = order[2 * j + 1];


        for (size_t jj = j + 1; jj < nchunk; jj++)
        {
            min_tmp = order[2 * jj];
            max_tmp = order[2 * jj + 1];


            if (((min >= min_tmp ) && (min <= max_tmp)) ||
                ((max >= min_tmp ) && (max <= max_tmp)))
            {
                sprintf(err_msg, "Chunks \"%lu, %lu\" and \"%lu, %lu\" "
                                 "overlap.", min, max, min_tmp, max_tmp);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG, err_msg);
                return -4;
            }
        }
    }


    /* Check that the chunks leave no gaps.  For instance, with "nmax = 10" and
     * "nchunk = 3", the "order" set to "0, 0, 2, 5, 6, 10" leaves a gap, as no
     * chunk contains the spherical harmonic order "1". */
    if (check_gaps)
    {
        _Bool match;
        unsigned long m = 0;
        for (size_t j = 0; j < nchunk; j++)
        {
            match = 0;
            for (size_t jj = 0; jj < nchunk; jj++)
            {
                if (order[2 * jj] == m)
                {
                    m     = order[2 * jj + 1] + 1;
                    match = 1;
                    break;
                }
            }


            if (!match)
            {
                sprintf(err_msg, "Couldn't find spherical harmonic order "
                                 "\"%lu\" in any chunk.", m);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG, err_msg);
                return -5;
            }
        }


        /* The value of "m" must now be "nmax + 1" */
        if (m != (nmax + 1))
        {
            sprintf(err_msg, "Couldn't find spherical harmonic order "
                             "\"nmax = %lu\" in any chunk.", nmax);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            return -6;
        }
    }


    return 0;
}
