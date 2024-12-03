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
#include "../src/prec.h"
#include "../src/err/err_is_null_ptr.h"
#include "../src/mpi/mpi_size_t.h"
#include "../src/mpi/mpi_allequal.h"
#include "parameters.h"
#include "error_messages.h"
#include "check_struct.h"
#include "partition_interval.h"
#include "shc_touch_array_elements.h"
#include "mpi_shc_distribute.h"
#include "check_mpi_shc.h"
#include "check_mpi_shc_alloc.h"
/* ------------------------------------------------------------------------- */






long int check_mpi_shc_alloc(CHARM(shc) *(*mpi_shc_alloc)(unsigned long,
                                                          REAL,
                                                          REAL,
                                                          size_t,
                                                         const unsigned long *,
                                                          MPI_Comm,
                                                          CHARM(err) *))
{
    /* --------------------------------------------------------------------- */
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    long int e = 0;
    REAL mu    = PREC(1.1);
    REAL r     = PREC(2.2);


    char func[NSTR_SHORT];
    if (mpi_shc_alloc == CHARM(mpi_shc_malloc))
        sprintf(func, "mpi_shc_malloc");
    else if (mpi_shc_alloc == CHARM(mpi_shc_calloc))
        sprintf(func, "mpi_shc_calloc");


    char func_call_str[NSTR_LONG];
    CHARM(shc) *shcs = NULL;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(mpi_err_init)(comm);
    if (CHARM(err_is_null_ptr)(err, 1, comm))
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    unsigned long nmax;
    size_t local_nchunk;
    for (int x = 0; x < 2; x++)
    {
        /* If "x == 0", all spherical harmonic coefficients will be distributed
         * on one MPI process only.  If "x == 1", the coefficients are
         * distributed across all MPI processes. */
    for (nmax = NMAX_MPI; nmax <= (NMAX_MPI + 10); nmax++)
    {
    for (local_nchunk = 1; local_nchunk <= NCHUNK_MAX; local_nchunk++)
    {
        unsigned long *local_order = mpi_shc_distribute(nmax, rank,
                                                        ((x == 0) &&
                                                         (rank == 0)) ?
                                                        1 : size,
                                                        local_nchunk);
        size_t nchunk_tmp = ((x == 0) && (rank != 0)) ? 0 : local_nchunk;


        shcs = mpi_shc_alloc(nmax, mu, r, nchunk_tmp, local_order, comm, err);
        CHARM(err_handler)(err, 0);


        sprintf(func_call_str, "%s(%lu, " FORMAT ", " FORMAT ", %zu, {",
                func, nmax, mu, r, local_nchunk);
        for (size_t j = 0; j < 2 * local_nchunk - 1; j++)
            sprintf(func_call_str + strlen(func_call_str), "%lu, ",
                    local_order[j]);
        sprintf(func_call_str + strlen(func_call_str), "%lu}, "
                "MPI_COMM_WORLD, err)", local_order[2 * local_nchunk - 1]);


        e += check_mpi_shc(shcs, func_call_str, err);


        e += check_struct_size_t(shcs->local_nchunk, nchunk_tmp, NEQ, VALID,
                              func_call_str,
                              "returned wrong value of \"local_nchunk\"");


        e += check_struct__Bool(shcs->owner, 1, NEQ, VALID, func_call_str,
                                "returned wrong value of \"owner\"");



        for (size_t j = 0; j < nchunk_tmp; j++)
        {
            e += check_struct_ulong(shcs->local_order[2 * j],
                                    local_order[2 * j], NEQ, VALID,
                                    func_call_str,
                                    "returned wrong value of an element of "
                                    "\"local_order\"");
            e += check_struct_ulong(shcs->local_order[2 * j + 1],
                                    local_order[2 * j + 1], NEQ, VALID,
                                    func_call_str,
                                    "returned wrong value of an element of "
                                    "\"local_order\"");
        }


        shc_touch_array_elements(shcs);
        CHARM(shc_free)(shcs);
        free(local_order);
    }
    }
    }
    /* --------------------------------------------------------------------- */


    /* Check that "mpi_shc_alloc" throws an error if "err" is not distributed
     * */
    /* --------------------------------------------------------------------- */
    {
    CHARM(err) *err_d0 = CHARM(err_init)();
    if (err_d0 == NULL)
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    local_nchunk = 1;
    unsigned long local_order[2] = {0, NMAX_MPI};
    shcs = mpi_shc_alloc(NMAX_MPI, mu, r, local_nchunk, local_order, comm,
                         err_d0);
    if (err_d0->code == CHARM_SUCCESS)
    {
        char err_msg[NSTR_LONG];
        sprintf(err_msg, "\"%s\" incorrectly accepted non-distributed "
                         "\"charm" CHARM_SUFFIX "_err\" structure", func);
        fprintf(stderr, err_msg);
        e += 1;
    }
    /* Call error handler only if "mpi_shc_alloc" wrongly did not report any
     * error */
    if (err_d0->code == CHARM_SUCCESS)
        CHARM(err_handler)(err_d0, 0);


    CHARM(err_free)(err_d0);
    CHARM(shc_free)(shcs);
    }
    /* --------------------------------------------------------------------- */


    /* Check that wrong values of "local_order" throw errors */
    /* --------------------------------------------------------------------- */
#   define NWRONG_CASES 7
#   define NCHUNK_TMP 3
    unsigned long local_order[NWRONG_CASES * 2 *
                              NCHUNK_TMP] = { 0, 30, 40, 50, 60, 101,
                                              0, 30, 40, 50, 60, 8,
                                             10, 30, 40, 50, 60, 100,
                                              0, 30, 30, 50, 60, 100,
                                              0, 40, 30, 50, 60, 100,
                                              0, 30,  0, 30, 40, 100,
                                             30,  0, 40, 50, 60, 100};


    for (int i = 0; i < NWRONG_CASES; i++)
    {
        shcs = mpi_shc_alloc(NMAX_MPI, mu, r, NCHUNK_TMP,
                             &local_order[i * 2 * NCHUNK_TMP], comm, err);


        sprintf(func_call_str, "%s(%lu, " FORMAT ", " FORMAT ", %zu, {",
                func, NMAX_MPI, mu, r, local_nchunk);
        for (size_t j = 0; j < 2 * local_nchunk - 1; j++)
            sprintf(func_call_str + strlen(func_call_str), "%lu, ",
                    local_order[j]);
        sprintf(func_call_str + strlen(func_call_str), "%lu}",
                local_order[2 * local_nchunk - 1]);


        /* Call error handler only if "mpi_shc_alloc" wrongly did not report
         * any error */
        if (err->code == CHARM_SUCCESS)
        {
            char err_msg[2 * NSTR_LONG];
            sprintf(err_msg, "\"%s\" incorrectly accepted the following "
                             "function call \"%s\".", func, func_call_str);
            fprintf(stderr, err_msg);
            e += 1;


            CHARM(err_handler)(err, 0);
        }


        CHARM(shc_free)(shcs);
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}

