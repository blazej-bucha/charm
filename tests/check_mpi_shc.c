/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdint.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../src/prec.h"
#include "../src/mpi/mpi_allequal.h"
#include "../src/mpi/mpi_size_t.h"
#include "../src/mpi/mpi_real.h"
#include "parameters.h"
#include "check_struct.h"
#include "check_mpi_shc.h"
/* ------------------------------------------------------------------------- */






long int check_mpi_shc(CHARM(shc) *shcs,
                       const char *func_call_str,
                       CHARM(err) *err)
{
    long int e = 0;


    e += check_struct_ptr(shcs, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer");


    long int etmp = 0;
    etmp = !CHARM(mpi_allequal_ulong)(shcs->nmax, shcs->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr,
                "\"shcs->nmax\" does not match across MPI processes\n");
    e += etmp;


    etmp = !CHARM(mpi_allequal_real)(shcs->mu, shcs->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr, "\"shcs->mu\" does not match across MPI processes\n");
    e += etmp;


    etmp = !CHARM(mpi_allequal_real)(shcs->r, shcs->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr, "\"shcs->r\" does not match across MPI processes\n");
    e += etmp;


    etmp = !CHARM(mpi_allequal_size_t)(shcs->nc, shcs->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr, "\"shcs->nc\" does not match across MPI processes\n");
    e += etmp;


    etmp = !CHARM(mpi_allequal_size_t)(shcs->ns, shcs->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr, "\"shcs->ns\" does not match across MPI processes\n");
    e += etmp;


    etmp = !CHARM(mpi_allequal__Bool)(shcs->owner, shcs->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr,
                "\"shcs->owner\" does not match across MPI processes\n");
    e += etmp;


    e += check_struct__Bool(shcs->distributed, 1, NEQ, VALID, func_call_str,
                            "returned wrong value of \"distributed\"");


    etmp = !CHARM(mpi_allequal__Bool)(shcs->distributed, shcs->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr,
                "\"shcs->distributed\" does not match across MPI processes\n");
    e += etmp;


    size_t local_nc_sum;
    MPI_Allreduce(&(shcs->local_nc), &local_nc_sum, 1, CHARM_MPI_SIZE_T,
                  MPI_SUM, shcs->comm);
    e += check_struct_size_t(shcs->nc, local_nc_sum, NEQ, VALID, func_call_str,
                             "returned wrong value of \"local_nc\"");


    size_t local_ns_sum;
    MPI_Allreduce(&(shcs->local_ns), &local_ns_sum, 1, CHARM_MPI_SIZE_T,
                  MPI_SUM, shcs->comm);
    e += check_struct_size_t(shcs->ns, local_ns_sum, NEQ, VALID, func_call_str,
                             "returned wrong value of \"local_ns\"");


    int mpi_compare;
    MPI_Comm_compare(MPI_COMM_WORLD, shcs->comm, &mpi_compare);
    e += check_struct_int(mpi_compare, MPI_IDENT, NEQ, VALID, func_call_str,
                          "returned wrong value of \"comm\"");


    return e;
}

