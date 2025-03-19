/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../src/prec.h"
#include "../src/crd/crd_point_isGLGrid.h"
#include "../src/crd/crd_point_isDHGrid.h"
#include "../src/err/err_is_null_ptr.h"
#include "../src/err/err_propagate.h"
#include "../src/mpi/mpi_size_t.h"
#include "../src/mpi/mpi_real.h"
#include "../src/mpi/mpi_count_aint.h"
#include "../src/mpi/mpi_gatherv.h"
#include "../src/mpi/mpi_size_t2charm_mpi_count.h"
#include "parameters.h"
#include "error_messages.h"
#include "cmp_vals.h"
#include "cmp_arrays.h"
#include "partition_interval.h"
#include "check_struct.h"
#include "mpi_crd_point_distribute.h"
#include "check_mpi_crd_point.h"
#include "check_mpi_crd_point_quad.h"
/* ------------------------------------------------------------------------- */






long int check_mpi_crd_point_quad(CHARM(point) *(*mpi_quad)(unsigned long,
                                                            REAL,
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
    REAL r     = PREC(1.1);


    char func[NSTR_SHORT];
    int grd_type = INT_MAX;
    CHARM(point) *(*quad)(unsigned long, REAL) = NULL;
    void (*shape)(unsigned long, size_t *, size_t *) = NULL;;
    if (mpi_quad == CHARM(mpi_crd_point_gl))
    {
        grd_type = CHARM_CRD_POINT_GRID_GL;
        quad  = CHARM(crd_point_gl);
        shape = CHARM(crd_point_gl_shape);
        snprintf(func, NSTR_SHORT, "mpi_crd_point_gl");
    }
    else if (mpi_quad == CHARM(mpi_crd_point_dh1))
    {
        grd_type = CHARM_CRD_POINT_GRID_DH1;
        quad  = CHARM(crd_point_dh1);
        shape = CHARM(crd_point_dh1_shape);
        snprintf(func, NSTR_SHORT, "mpi_crd_point_dh1");
    }
    else if (mpi_quad == CHARM(mpi_crd_point_dh2))
    {
        grd_type = CHARM_CRD_POINT_GRID_DH2;
        quad  = CHARM(crd_point_dh2);
        shape = CHARM(crd_point_dh2_shape);
        snprintf(func, NSTR_SHORT, "mpi_crd_point_dh2");
    }


    char func_call_str[NSTR_LONG];
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
    CHARM(point) *grd_d0 = NULL;
    CHARM(point) *grd_d1 = NULL;
    size_t nlat, nlon;


    for (unsigned long nmax = 0; nmax <= NMAX_MPI; nmax++)
    {
        grd_d0 = quad(nmax, r);
        if (grd_d0 == NULL)
        {
            fprintf(stderr, ERR_MSG_POINT);
            exit(CHARM_FAILURE);
        }


        /* Compute the distributed quadrature grid */
        /* ................................................................. */
        shape(nmax, &nlat, &nlon);


        size_t local_nlat, local_0_start;
        mpi_crd_point_distribute(nlat, grd_d0->type, rank, size,
                                 &local_nlat, &local_0_start);


        grd_d1 = mpi_quad(nmax, r, local_nlat, local_0_start, comm, err);
        CHARM(err_handler)(err, 0);


        snprintf(func_call_str, NSTR_LONG,
                                "%s(%lu, " REAL_PRINT_FORMAT
                                ", %zu, %zu, MPI_COMM_WORLD, err)",
                                func, nmax, r, local_nlat, local_0_start);


        e += check_struct_ptr(grd_d1, NULL, EQ, VALID, func_call_str,
                              "returned NULL pointer");
        /* ................................................................. */


        /* Gather data from "grd_d1" to compare it with "grd_d0" */
        /* ................................................................. */
        REAL *lat = (REAL *)malloc(grd_d0->nlat * sizeof(REAL));
        if (lat == NULL)
        {
            fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
            exit(CHARM_FAILURE);
        }


        REAL *lon = (REAL *)malloc(grd_d0->nlon * sizeof(REAL));
        if (lon == NULL)
        {
            fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
            exit(CHARM_FAILURE);
        }


        REAL *r = (REAL *)malloc(grd_d0->nlat * sizeof(REAL));
        if (r == NULL)
        {
            fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
            exit(CHARM_FAILURE);
        }


        REAL *w = (REAL *)malloc(grd_d0->nlat * sizeof(REAL));
        if (w == NULL)
        {
            fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
            exit(CHARM_FAILURE);
        }


        CHARM_MPI_COUNT *recvcounts =
                     (CHARM_MPI_COUNT *)malloc(size * sizeof(CHARM_MPI_COUNT));
        if (recvcounts == NULL)
        {
            fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
            exit(CHARM_FAILURE);
        }
        CHARM_MPI_COUNT local_nlat_c =
                         CHARM(mpi_size_t2charm_mpi_count)(grd_d1->local_nlat,
                                                           err);
        CHARM(err_handler)(err, 1);


        MPI_Allgather(&local_nlat_c, 1, CHARM_MPI_COUNT_DATATYPE,
                      recvcounts, 1, CHARM_MPI_COUNT_DATATYPE, comm);


        CHARM_MPI_AINT *displs =
                       (CHARM_MPI_AINT *)malloc(size * sizeof(CHARM_MPI_AINT));
        if (displs == NULL)
        {
            fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
            exit(CHARM_FAILURE);
        }
        displs[0] = 0;
        for (int i = 1; i < size; i++)
            displs[i] = displs[i - 1] + recvcounts[i - 1];


        CHARM_MPI_GATHERV(grd_d1->lat, grd_d1->local_nlat, CHARM_MPI_REAL,
                          lat, recvcounts, displs, CHARM_MPI_REAL, 0, comm);
        CHARM_MPI_GATHERV(grd_d1->r, grd_d1->local_nlat, CHARM_MPI_REAL,
                          r, recvcounts, displs, CHARM_MPI_REAL, 0, comm);
        CHARM_MPI_GATHERV(grd_d1->w, grd_d1->local_nlat, CHARM_MPI_REAL,
                          w, recvcounts, displs, CHARM_MPI_REAL, 0, comm);
        /* ................................................................. */


        /* ................................................................. */
        size_t x = (size_t)CHARM(crd_point_isDHGrid)(grd_d1->type);
        size_t global_north, global_south, local_south;
        for (size_t i = 0; i < (grd_d1->local_nlat + 1) / 2; i++)
        {
            global_north = grd_d1->local_0_start + i;
            e += cmp_vals_real(grd_d1->lat[i], grd_d0->lat[global_north],
                               CHARM(glob_threshold));


            e += cmp_vals_real(grd_d1->r[i], grd_d0->r[global_north],
                               CHARM(glob_threshold));


            e += cmp_vals_real(grd_d1->w[i], grd_d0->w[global_north],
                               CHARM(glob_threshold));


            /* Do not use the equatorial symmetry at the north pole for
             * Driscoll--Healy grids */
            if (x && (grd_d1->local_0_start == 0) && (i == 0))
                continue;


            size_t p = 0;
            /* If the local portion of the Driscoll--Healy grid contains the
             * north-pole, the symmetric latitudes on the south hemisphere have
             * to be switched by "1", as there is no south pole */
            if (x && (grd_d1->local_0_start == 0))
                p = 1;


            local_south = grd_d1->local_nlat - 1 - i + p;
            global_south = grd_d0->nlat - 1 - (grd_d1->local_0_start + i - x);
            e += cmp_vals_real(grd_d1->lat[local_south],
                               grd_d0->lat[global_south],
                               CHARM(glob_threshold));


            e += cmp_vals_real(grd_d1->r[local_south], grd_d0->r[global_south],
                               CHARM(glob_threshold));


            e += cmp_vals_real(grd_d1->w[local_south], grd_d0->w[global_south],
                               CHARM(glob_threshold));
        }
        e += cmp_arrays(grd_d1->lon, grd_d0->lon, grd_d0->nlon,
                        CHARM(glob_threshold));
        /* ................................................................. */


        /* ................................................................. */
        e += check_mpi_crd_point(grd_d1, func_call_str, err);


        e += check_struct__Bool(grd_d1->owner, 1, NEQ, VALID, func_call_str,
                                "returned wrong value of \"owner\"");
        /* ................................................................. */


        /* ................................................................. */
        CHARM(crd_point_free)(grd_d0);
        CHARM(crd_point_free)(grd_d1);
        free(lat);
        free(lon);
        free(r);
        free(w);
        free(recvcounts);
        free(displs);
        /* ................................................................. */
    }
    /* --------------------------------------------------------------------- */


    /* Check that wrong total number of points returns an error */
    /* --------------------------------------------------------------------- */
    {
    shape(NMAX_MPI, &nlat, &nlon);


    size_t local_nlat, local_0_start;
    mpi_crd_point_distribute(nlat, grd_type, rank, size,
                             &local_nlat, &local_0_start);


    if (rank == 0)
        local_nlat += 1;


    grd_d1 = mpi_quad(NMAX_MPI, r, local_nlat, local_0_start, comm, err);


    snprintf(func_call_str, NSTR_LONG,
                            "%s(%lu, " REAL_PRINT_FORMAT
                            ", %zu, %zu, MPI_COMM_WORLD, err)",
                            func, NMAX_MPI, r, local_nlat, local_0_start);


    e += check_struct_ptr(grd_d1, NULL, NEQ, INVALID, func_call_str,
                          "didn't return NULL pointer");
    if (err->code == CHARM_SUCCESS)
        CHARM(err_handler)(err, 0);


    CHARM(crd_point_free)(grd_d1);
    }
    /* --------------------------------------------------------------------- */


    /* Check some of the wrong distributions of latitudinal chunks that have to
     * return NULL pointer. */
    /* --------------------------------------------------------------------- */
    if (size > 1)
    {
        shape(NMAX_MPI, &nlat, &nlon);


        size_t local_nlat, local_0_start;
        mpi_crd_point_distribute(nlat, grd_type, rank, size,
                                 &local_nlat, &local_0_start);


        if (CHARM(crd_point_isGLGrid)(grd_type))
        {
            if (local_nlat % 2)
                local_nlat -= 1;
            if (rank == 0)
                local_nlat += 1;

        }
        else if (CHARM(crd_point_isDHGrid)(grd_type))
        {
            if (rank == 0)
                local_nlat += 1;
            else if (rank == (size - 1))
                local_nlat -= 1;
        }


        grd_d1 = mpi_quad(NMAX_MPI, r, local_nlat, local_0_start, comm, err);


        snprintf(func_call_str, NSTR_LONG,
                                "%s(%lu, " REAL_PRINT_FORMAT
                                ", %zu, %zu, MPI_COMM_WORLD, err)",
                                func, NMAX_MPI, r, local_nlat, local_0_start);


        e += check_struct_ptr(grd_d1, NULL, NEQ, INVALID, func_call_str,
                              "didn't return NULL pointer");
        if (err->code == CHARM_SUCCESS)
            CHARM(err_handler)(err, 0);


        CHARM(crd_point_free)(grd_d1);
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


    shape(NMAX_MPI, &nlat, &nlon);


    size_t local_nlat, local_0_start;
    mpi_crd_point_distribute(nlat, grd_type, rank, size,
                             &local_nlat, &local_0_start);


    grd_d1 = mpi_quad(NMAX_MPI, r, local_nlat, local_0_start, comm, err_d0);


    if (err_d0->code == CHARM_SUCCESS)
    {
        char err_msg[NSTR_LONG];
        snprintf(err_msg, NSTR_LONG,
                          "\"%s\" incorrectly accepted non-distributed "
                          "\"charm" CHARM_SUFFIX "_err\" structure", func);
        fprintf(stderr, err_msg);
        e += 1;


        CHARM(err_handler)(err_d0, 0);
    }


    CHARM(crd_point_free)(grd_d1);
    CHARM(err_free)(err_d0);
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    return e;
}

