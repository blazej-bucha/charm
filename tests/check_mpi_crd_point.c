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
#include "../src/crd/crd_point_isQuadGrid.h"
#include "parameters.h"
#include "check_struct.h"
#include "point_touch_array_elements.h"
#include "check_mpi_crd_point.h"
/* ------------------------------------------------------------------------- */






long int check_mpi_crd_point(CHARM(point) *pnt,
                             const char *func_call_str,
                             CHARM(err) *err)
{
    long int e = 0;


    e += check_struct_ptr(pnt, NULL, EQ, VALID, func_call_str,
                          "returned NULL pointer");


    e += check_struct__Bool(pnt->distributed, 1, NEQ, VALID, func_call_str,
                            "returned wrong value of distributed\"");


    long int etmp = 0;
    etmp = !CHARM(mpi_allequal_int)(pnt->type, pnt->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr, "\"pnt->type\" does not match across MPI processes\n");
    e += etmp;


    etmp = !CHARM(mpi_allequal_size_t)(pnt->nlat, pnt->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr, "\"pnt->nlat\" does not match across MPI processes\n");
    e += etmp;


    etmp = !CHARM(mpi_allequal_size_t)(pnt->nlon, pnt->comm, err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr, "\"pnt->nlon\" does not match across MPI processes\n");
    e += etmp;


    etmp = !CHARM(mpi_allequal_size_t)(pnt->npoint, pnt->comm,
                                       err);
    CHARM(err_handler)(err, 1);
    if (etmp)
        fprintf(stderr,
                "\"pnt->npoint\" does not match across MPI processes\n");
    e += etmp;


    size_t tmp;
    MPI_Allreduce(&(pnt->local_nlat), &tmp, 1, CHARM_MPI_SIZE_T, MPI_SUM,
                  pnt->comm);
    e += check_struct_size_t(pnt->nlat, tmp, NEQ, VALID, func_call_str,
                             "returned wrong value of \"local_nlat\"");


    MPI_Allreduce(&(pnt->local_npoint), &tmp, 1, CHARM_MPI_SIZE_T, MPI_SUM,
                  pnt->comm);
    e += check_struct_size_t(pnt->npoint, tmp, NEQ, VALID, func_call_str,
                             "returned wrong value of \"local_npoint\"");


    int mpi_compare;
    MPI_Comm_compare(MPI_COMM_WORLD, pnt->comm, &mpi_compare);
    e += check_struct_int(mpi_compare, MPI_IDENT, NEQ, VALID, func_call_str,
                          "returned wrong value of \"comm\"");


    point_touch_array_elements(pnt);


    return e;
}

