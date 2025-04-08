/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/err/err_is_null_ptr.h"
#include "parameters.h"
#include "error_messages.h"
#include "cmp_vals.h"
#include "check_mpi_shc_local_ncs.h"
/* ------------------------------------------------------------------------- */






long int check_mpi_shc_local_ncs(void)
{
    /* We use the same code to check "mpi_shc_local_ncs" on all processes.
     * This is because an automatic split of coefficients based on "size" and
     * "rank" would mean we cannot hard code the reference values. */
    /* --------------------------------------------------------------------- */
    MPI_Comm comm = MPI_COMM_WORLD;
    CHARM(err) *err = CHARM(mpi_err_init)(comm);
    if (CHARM(err_is_null_ptr)(err, 1, comm))
    {
        fprintf(stderr, "%s", ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    long int e = 0;
    unsigned long nmax = 10;


    /* Zero local chunks means zero local coefficients */
    {
    unsigned long local_order[2] = {0, nmax};
    size_t local_ncs = CHARM(mpi_shc_local_ncs)(nmax, 0, local_order, err);
    CHARM(err_handler)(err, 0);
    e += cmp_vals_ulong(local_ncs, 0);
    }


    /* One full chunk from "0" to "nmax" */
    {
    unsigned long local_order[2] = {0, nmax};
    size_t local_ncs = CHARM(mpi_shc_local_ncs)(nmax, 1, local_order, err);
    CHARM(err_handler)(err, 0);
    e += cmp_vals_ulong(local_ncs, 66);
    }


    /* One partial chunk */
    {
    unsigned long local_order[2] = {0, 0};
    size_t local_ncs = CHARM(mpi_shc_local_ncs)(nmax, 1, local_order, err);
    CHARM(err_handler)(err, 0);
    e += cmp_vals_ulong(local_ncs, 11);
    }
    {
    unsigned long local_order[2] = {2, 4};
    size_t local_ncs = CHARM(mpi_shc_local_ncs)(nmax, 1, local_order, err);
    CHARM(err_handler)(err, 0);
    e += cmp_vals_ulong(local_ncs, 24);
    }
    {
    unsigned long local_order[2] = {nmax, nmax};
    size_t local_ncs = CHARM(mpi_shc_local_ncs)(nmax, 1, local_order, err);
    CHARM(err_handler)(err, 0);
    e += cmp_vals_ulong(local_ncs, 1);
    }


    /* Two chunks */
    {
    unsigned long local_order[4] = {0, 3, 4, nmax};
    size_t local_ncs = CHARM(mpi_shc_local_ncs)(nmax, 2, local_order, err);
    CHARM(err_handler)(err, 0);
    e += cmp_vals_ulong(local_ncs, 66);
    }


    /* Two partial chunks */
    {
    unsigned long local_order[4] = {1, 2, 3, 3};
    size_t local_ncs = CHARM(mpi_shc_local_ncs)(nmax, 2, local_order, err);
    CHARM(err_handler)(err, 0);
    e += cmp_vals_ulong(local_ncs, 27);
    }


    /* Three chunks */
    {
    unsigned long local_order[6] = {0, 3, 4, 6, 7, nmax};
    size_t local_ncs = CHARM(mpi_shc_local_ncs)(nmax, 3, local_order, err);
    CHARM(err_handler)(err, 0);
    e += cmp_vals_ulong(local_ncs, 66);
    }


    /* Three partial chunks */
    {
    unsigned long local_order[6] = {1, 1, 2, 2, 3, 3};
    size_t local_ncs = CHARM(mpi_shc_local_ncs)(nmax, 3, local_order, err);
    CHARM(err_handler)(err, 0);
    e += cmp_vals_ulong(local_ncs, 27);
    }
    /* --------------------------------------------------------------------- */


    /* Check that "mpi_shc_local_ncs" throws an error if "err" is not
     * distributed */
    /* --------------------------------------------------------------------- */
    {
    CHARM(err) *err_d0 = CHARM(err_init)();
    if (err_d0 == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    nmax = 100;
    unsigned long local_order[2] = {0, nmax};
    CHARM(mpi_shc_local_ncs)(nmax, 1, local_order, err_d0);
    if (err_d0->code == CHARM_SUCCESS)
    {
        fprintf(stderr, "%s",
                "charm" CHARM_SUFFIX "_mpi_local_shcs_local_ncs "
                "incorrectly accepted non-distributed "
                "\"charm" CHARM_SUFFIX "_err\" structure");
        e += 1;
    }
    /* Call error handler only if "mpi_shc_alloc" wrongly did not report any
     * error */
    if (err_d0 == CHARM_SUCCESS)
        CHARM(err_handler)(err_d0, 0);


    CHARM(err_free)(err_d0);
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}

