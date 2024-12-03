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
#include "parameters.h"
#include "error_messages.h"
#include "check_struct.h"
#include "check_mpi_crd_point.h"
#include "check_mpi_crd_point_init.h"
/* ------------------------------------------------------------------------- */






/* Checks "mpi_crd_point_init".  Assumes that "mpi_crd_point_calloc" and
 * "mpi_crd_point_malloc" have already been tested, as it checks only the
 * features specifically related to "mpi_shc_init". */
long int check_mpi_crd_point_init(void)
{
    /* --------------------------------------------------------------------- */
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    long int e = 0;
    const size_t local_nlat = 5;
    const size_t nlon = 5;


    const char func[NSTR_SHORT] = "mpi_crd_point_init";
    int type = CHARM_CRD_POINT_SCATTERED;
    char func_call_str[NSTR_LONG];
    CHARM(point) *pnt = NULL;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(mpi_err_init)(comm);
    if (CHARM(err_is_null_ptr)(err, 1, comm))
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    size_t local_0_start = rank * local_nlat;


    REAL *lat = (REAL *)malloc(local_nlat * sizeof(REAL));
    if (lat == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    REAL *lon = (REAL *)malloc(nlon * sizeof(REAL));
    if (lon == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    REAL *r = (REAL *)malloc(local_nlat * sizeof(REAL));
    if (r == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }


    pnt = CHARM(mpi_crd_point_init)(type, local_nlat, nlon, local_0_start,
                                    lat, lon, r, comm, err);
    CHARM(err_handler)(err, 0);


    sprintf(func_call_str, "%s(%d, %zu, %zu, %zu, MPI_COMM_WORLD, "
                           "err)", func, type, local_nlat, nlon,
                           local_0_start);


    e += check_mpi_crd_point(pnt, func_call_str, err);


    e += check_struct__Bool(pnt->owner, 0, NEQ, VALID, func_call_str,
                            "returned wrong value of \"owner\"");


    e += check_struct_ptr(pnt->lat, lat, NEQ, VALID, func_call_str,
                          "returned wrong value of \"lat\"");


    e += check_struct_ptr(pnt->lon, lon, NEQ, VALID, func_call_str,
                          "returned wrong value of \"lon\"");


    e += check_struct_ptr(pnt->r, r, NEQ, VALID, func_call_str,
                          "returned wrong value of \"r\"");


    CHARM(crd_point_free)(pnt);
    /* --------------------------------------------------------------------- */


    /* Check that an error is thrown if "err" is not distributed */
    /* --------------------------------------------------------------------- */
    {
    CHARM(err) *err_d0 = CHARM(err_init)();
    if (err_d0 == NULL)
    {
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    pnt = CHARM(mpi_crd_point_init)(type, 1, 1, rank, lat, lon, r, comm,
                                    err_d0);


    if (err_d0->code == CHARM_SUCCESS)
    {
        char err_msg[NSTR_LONG];
        sprintf(err_msg, "\"%s\" incorrectly accepted non-distributed "
                         "\"charm" CHARM_SUFFIX "_err\" structure", func);
        fprintf(stderr, err_msg);
        e += 1;


        CHARM(err_handler)(err_d0, 0);
    }


    CHARM(crd_point_free)(pnt);
    CHARM(err_free)(err_d0);
    }
    /* --------------------------------------------------------------------- */



    /* --------------------------------------------------------------------- */
    free(lat);
    free(lon);
    free(r);
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}

