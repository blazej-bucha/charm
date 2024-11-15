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
#include "check_mpi_crd_point_alloc.h"
/* ------------------------------------------------------------------------- */






/* Number of point types to be checked */
/* ------------------------------------------------------------------------- */
#undef POINT_TYPES
#define POINT_TYPES (2)
/* ------------------------------------------------------------------------- */






long int check_mpi_crd_point_alloc(CHARM(point) *(*mpi_alloc)(int,
                                                              size_t,
                                                              size_t,
                                                              size_t,
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
    const size_t local_nlat = 5;
    const size_t nlon = 10;


    char func[NSTR_SHORT];
    if (mpi_alloc == CHARM(mpi_crd_point_malloc))
        sprintf(func, "mpi_crd_point_malloc");
    else if (mpi_alloc == CHARM(mpi_crd_point_calloc))
        sprintf(func, "mpi_crd_point_calloc");


    int type;
    /* We do not check here the GL, DH1 and DH2 point types, as they require
     * a different treatment of "local_0_start".  The check for GL, DH1 and DH2
     * point types are performed by "check_mpi_crd_point_quad.c" */
    int types[POINT_TYPES] = {CHARM_CRD_POINT_SCATTERED,
                              CHARM_CRD_POINT_GRID};


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
    for (int t = 0; t < POINT_TYPES; t++)
    {
        type = types[t];


        for (size_t i = 0; i < local_nlat; i++)
        {
            size_t local_0_start = rank * i;


            for (size_t j = 0; j < nlon; j++)
            {
                pnt = mpi_alloc(type, i, j, local_0_start, comm, err);
                if ((type == CHARM_CRD_POINT_SCATTERED) && (i != j))
                    CHARM(err_reset)(err);
                else
                    CHARM(err_handler)(err, 0);


                sprintf(func_call_str, "%s(%d, %zu, %zu, %zu, MPI_COMM_WORLD, "
                                       "err)", func, type, i, j,
                                       local_0_start);


                if ((type == CHARM_CRD_POINT_SCATTERED) && (i != j))
                    e += check_struct_ptr(pnt, NULL, NEQ, INVALID,
                                          func_call_str,
                                          "didn't returned NULL pointer");
                else
                    e += check_struct_ptr(pnt, NULL, EQ, VALID, func_call_str,
                                          "returned NULL pointer");


                if (pnt != NULL)
                {
                    e += check_mpi_crd_point(pnt, func_call_str, err);


                    e += check_struct__Bool(pnt->owner, 1, NEQ, VALID,
                                            func_call_str,
                                            "returned wrong value of "
                                            "\"owner\"");
                }


                CHARM(crd_point_free)(pnt);
            }
        }
    }
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


    pnt = mpi_alloc(CHARM_CRD_POINT_SCATTERED, 1, 1, rank, comm, err_d0);


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
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}

