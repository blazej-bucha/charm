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
#include "partition_interval.h"
#include "shc_touch_array_elements.h"
#include "check_struct.h"
#include "check_mpi_shc.h"
#include "check_mpi_shc_init.h"
/* ------------------------------------------------------------------------- */






/* Checks "mpi_shc_init".  Assumes that "mpi_shc_calloc", "mpi_shc_malloc"
 * "mpi_shc_local_ncs" have already been tested, as it checks only the features
 * specifically related to "mpi_shc_init". */
long int check_mpi_shc_init(void)
{
    /* --------------------------------------------------------------------- */
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    long int e                   = 0;
    REAL mu                      = PREC(1.1);
    REAL r                       = PREC(2.2);
    size_t local_nchunk          = 1;


    unsigned long start, end;
    partition_interval_ulong(0, NMAX_MPI, size, rank, &start, &end);
    unsigned long local_order[2] = {start, end};


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


    size_t local_ncs = CHARM(mpi_shc_local_ncs)(NMAX_MPI, local_nchunk,
                                                local_order, err);
    CHARM(err_handler)(err, 1);


    REAL *c = (REAL *)malloc(local_ncs * sizeof(REAL));
    if (c == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    REAL *s = (REAL *)malloc(local_ncs * sizeof(REAL));
    if (s == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    shcs = CHARM(mpi_shc_init)(NMAX_MPI, mu, r, c, s, local_nchunk,
                               local_order, comm, err);
    CHARM(err_handler)(err, 0);


    sprintf(func_call_str, "%s(%lu, " REAL_PRINT_FORMAT ", " REAL_PRINT_FORMAT
                           ", c, s, %zu, {",
            "charm" CHARM_SUFFIX "_mpi_shc_init", NMAX_MPI, mu, r, local_nchunk);
    for (size_t j = 0; j < 2 * local_nchunk - 1; j++)
        sprintf(func_call_str + strlen(func_call_str), "%lu, ",
                local_order[j]);
    sprintf(func_call_str + strlen(func_call_str), "%lu}, "
            "MPI_COMM_WORLD, err)", local_order[2 * local_nchunk - 1]);


    e += check_mpi_shc(shcs, func_call_str, err);


    e += check_struct_ptr(shcs, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer");


    e += check_struct__Bool(shcs->owner, 0, NEQ, VALID,
                            func_call_str,
                            "returned wrong value of \"owner\"");


    e += check_struct_ptr(shcs->c[shcs->local_order[0]], c, NEQ, VALID,
                          func_call_str, "returned wrong value of \"c\"");


    e += check_struct_ptr(shcs->s[shcs->local_order[0]], s, NEQ, VALID,
                          func_call_str, "returned wrong value of \"s\"");


    shc_touch_array_elements(shcs);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    free(c);
    free(s);
    CHARM(shc_free)(shcs);
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}

