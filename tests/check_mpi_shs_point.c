/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../src/prec.h"
#include "../src/err/err_is_null_ptr.h"
#include "../src/crd/crd_point_isSctr.h"
#include "../src/crd/crd_point_isCustGrid.h"
#include "../src/crd/crd_point_isQuadGrid.h"
#include "../src/crd/crd_point_isGLGrid.h"
#include "../src/crd/crd_point_isDH1Grid.h"
#include "../src/crd/crd_point_isDH2Grid.h"
#include "parameters.h"
#include "error_messages.h"
#include "partition_interval.h"
#include "generate_point.h"
#include "cmp_arrays.h"
#include "mpi_shc_distribute.h"
#include "mpi_crd_point_distribute.h"
#include "check_mpi_shs_point.h"
/* ------------------------------------------------------------------------- */






/* ------------------------------------------------------------------------- */
#undef NPOINT_TYPE
#define NPOINT_TYPE (5)
/* ------------------------------------------------------------------------- */






long int check_mpi_shs_point(void)
{
    /* --------------------------------------------------------------------- */
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    long int e = 0;
    const size_t local_npoint = 20;


    const REAL mu = PREC(1.1);
    const REAL r  = PREC(1.2);
    const REAL dr = PREC(0.1);


    CHARM(point) *pnt_d0 = NULL;
    CHARM(point) *pnt_d1 = NULL;


    CHARM(shc) *shcs_d0 = NULL;
    CHARM(shc) *shcs_d1 = NULL;


    int point_types_d1[NPOINT_TYPE] = {CHARM_CRD_POINT_SCATTERED,
                                       CHARM_CRD_POINT_GRID,
                                       CHARM_CRD_POINT_GRID_GL,
                                       CHARM_CRD_POINT_GRID_DH1,
                                       CHARM_CRD_POINT_GRID_DH2};


    /* Because of how we synthesize the reference values, the point types
     * corresponding to distributed quadrature grids are "CHARM_CRD_POINT_GRID"
     * */
    int point_types_d0[NPOINT_TYPE] = {CHARM_CRD_POINT_SCATTERED,
                                       CHARM_CRD_POINT_GRID,
                                       CHARM_CRD_POINT_GRID,
                                       CHARM_CRD_POINT_GRID,
                                       CHARM_CRD_POINT_GRID};
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(err) *err_d1 = CHARM(mpi_err_init)(comm);
    if (CHARM(err_is_null_ptr)(err_d1, 1, comm))
    {
        fprintf(stderr, "%s", ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    CHARM(err) *err_d0 = CHARM(err_init)();
    if (err_d0 == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */



    unsigned block_s = CHARM(glob_shs_block_lat_multiplier);


    for (int t = 0; t < NPOINT_TYPE; t++)
    {
    int type_d0 = point_types_d0[t];
    int type_d1 = point_types_d1[t];


    for (unsigned b = 0; b < 3; b++)
    {
    CHARM(glob_shs_block_lat_multiplier) = block_s + b;


    for (int x = 0; x < 2; x++)
    {
    /* If "x == 0", all spherical harmonic coefficients will be distributed on
     * one MPI process only.  If "x == 1", the coefficients are distributed
     * across all MPI processes. */


    for (int s = 0; s < 2; s++)
    {
    /* If "s == 0", the grids are symmetric with respect to the equator.  If "s
     * == 1", the grids are non-symmetric. */


    for (unsigned long nmax = NMAX_MPI; nmax <= NMAX_MPI + 2; nmax++)
    {
    for (size_t local_nchunk = 1; local_nchunk <= NCHUNK_MAX; local_nchunk++)
    {
    /* Non-distributed coefficients */
    /* --------------------------------------------------------------------- */
    shcs_d0 = CHARM(shc_calloc)(nmax, mu, r);
    if (shcs_d0 == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }


    /* Generate some fake coefficients.  We have already set all coefficients
     * to zero, so by setting "C00" to "1.0", we get a constant function "1.0"
     * over the entire sphere, to which the coefficients refer to.  Otherwise,
     * one has to be very careful to ensure that the fake coefficients generate
     * such a function, which does not vary too much with latitude and
     * longitude.  If the variations are larger than, say, 10 orders of
     * magnitude, this will not allow us to validate the synthesis with
     * sufficient accuracy. */
    shcs_d0->c[0][0] = PREC(1.0);
    /* --------------------------------------------------------------------- */


    /* Distributed coefficients */
    /* --------------------------------------------------------------------- */
    unsigned long *local_order = mpi_shc_distribute(nmax, rank,
                                                    ((x == 0) &&
                                                     (rank == 0)) ? 1 : size,
                                                    local_nchunk);
    size_t nchunk_tmp = ((x == 0) && (rank != 0)) ? 0 : local_nchunk;


    shcs_d1 = CHARM(mpi_shc_init)(nmax, mu, r, shcs_d0->c[local_order[0]],
                                  shcs_d0->s[local_order[0]], nchunk_tmp,
                                  local_order, comm, err_d1);
    CHARM(err_handler)(err_d1, 1);


    free(local_order);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    for (size_t i = 0; i < local_npoint; i += 3)
    {
        for (size_t j = i; j < local_npoint; j += 3)
        {
            if (CHARM(crd_point_isSctr)(type_d1) && (i != j))
                continue;
            else if (CHARM(crd_point_isQuadGrid)(type_d1) && (i != j))
                continue;


            if (CHARM(crd_point_isSctr)(type_d1) ||
                CHARM(crd_point_isCustGrid)(type_d1))
            {
                /* Symmetric and non-symmetric grids do not make sense with
                 * scattered points */
                if (CHARM(crd_point_isSctr)(type_d1) && (s == 1))
                    continue;


                pnt_d1 = CHARM(mpi_crd_point_malloc)(type_d1, i, j, i * rank,
                                                     comm, err_d1);
                CHARM(err_handler)(err_d1, 1);


                /* Now generate some latitudes and longitudes that will be
                 * different for each MPI process.  The spherical radius
                 * remains constant */
                REAL latstart = PI * ((REAL)(rank + 1) / (REAL)(size));
                REAL lonstart = (PREC(2.0) * PI) * ((REAL)(rank + 1) /
                                                    (REAL)(size));
                CHARM(generate_point)(pnt_d1, r + dr, latstart, lonstart);
                if ((s == 1) && (pnt_d1->local_nlat > 0))
                    pnt_d1->lat[0] -= (REAL)(BREAK_SYMM);
            }
            else if (CHARM(crd_point_isQuadGrid)(type_d1))
            {
                /* Quadrature grids are always symmetric and global */
                if (s == 1)
                    continue;


                size_t nlat, nlon, local_nlat, local_0_start;
                unsigned long nmax2 = nmax + i;  /* This does not make much
                                                  * sense, except it is
                                                  * a validation. */
                if (CHARM(crd_point_isGLGrid)(type_d1))
                {
                    CHARM(crd_point_gl_shape)(nmax2, &nlat, &nlon);
                    mpi_crd_point_distribute(nlat, type_d1, rank, size,
                                             &local_nlat, &local_0_start);
                    pnt_d1 = CHARM(mpi_crd_point_gl)(nmax2, r + dr, local_nlat,
                                                     local_0_start, comm,
                                                     err_d1);
                }
                else if (CHARM(crd_point_isDH1Grid)(type_d1))
                {
                    CHARM(crd_point_dh1_shape)(nmax2, &nlat, &nlon);
                    mpi_crd_point_distribute(nlat, type_d1, rank, size,
                                             &local_nlat, &local_0_start);
                    pnt_d1 = CHARM(mpi_crd_point_dh1)(nmax2, r + dr,
                                                      local_nlat,
                                                      local_0_start, comm,
                                                      err_d1);
                }
                else if (CHARM(crd_point_isDH2Grid)(type_d1))
                {
                    CHARM(crd_point_dh2_shape)(nmax2, &nlat, &nlon);
                    mpi_crd_point_distribute(nlat, type_d1, rank, size,
                                             &local_nlat, &local_0_start);
                    pnt_d1 = CHARM(mpi_crd_point_dh2)(nmax2, r + dr,
                                                      local_nlat,
                                                      local_0_start, comm,
                                                      err_d1);
                }
                CHARM(err_handler)(err_d1, 1);
            }


            pnt_d0 = CHARM(crd_point_init)(type_d0,
                                           pnt_d1->local_nlat,
                                           pnt_d1->local_nlon,
                                           pnt_d1->lat,
                                           pnt_d1->lon,
                                           pnt_d1->r);
            if (pnt_d0 == NULL)
            {
                fprintf(stderr, "%s", ERR_MSG_POINT);
                exit(CHARM_FAILURE);
            }


            REAL *f_d1 = (REAL *)malloc(pnt_d1->local_npoint * sizeof(REAL));
            if (f_d1 == NULL)
            {
                fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
                exit(CHARM_FAILURE);
            }


            REAL *f_d0 = (REAL *)malloc(pnt_d0->npoint * sizeof(REAL));
            if (f_d0 == NULL)
            {
                fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
                exit(CHARM_FAILURE);
            }


            CHARM(shs_point)(pnt_d0, shcs_d0, shcs_d0->nmax, f_d0, err_d0);
            CHARM(err_handler)(err_d0, 1);


            CHARM(shs_point)(pnt_d1, shcs_d1, shcs_d1->nmax, f_d1, err_d1);
            CHARM(err_handler)(err_d1, 1);


            /* Before calling "cmp_arrays", keep in mind that
             * "pnt_d1->local_npoint == pnt_d0->npoint".  Just to be sure, we
             * check this. */
            if (pnt_d1->local_npoint != pnt_d0->npoint)
            {
                fprintf(stderr, "%s",
                        "The number of local points of a distributed "
                        "point structure does not match the number of "
                        "points of a non-distributed structure.\n");
                exit(CHARM_FAILURE);
            }


            e += cmp_arrays(f_d0, f_d1, pnt_d1->local_npoint,
                            CHARM(glob_threshold2));


            CHARM(crd_point_free)(pnt_d0);
            CHARM(crd_point_free)(pnt_d1);
            free(f_d0);
            free(f_d1);
        }
    }


    CHARM(shc_free)(shcs_d0);
    CHARM(shc_free)(shcs_d1);
    }
    }
    }
    }
    }
    }
    /* --------------------------------------------------------------------- */


    /* Check that "charm_pnt", "charm_shc" and "charm_err" are distributed */
    /* --------------------------------------------------------------------- */
    /* "charm_point" is not distributed */
    /* ..................................................................... */
    {
    pnt_d0  = CHARM(crd_point_gl)(NMAX_MPI, r);
    if (pnt_d0 == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_POINT);
        exit(CHARM_FAILURE);
    }
    unsigned long local_order[2] = {0, NMAX_MPI};
    shcs_d1 = CHARM(mpi_shc_malloc)(NMAX_MPI, mu, r, 1 * (rank == 0),
                                    local_order, comm, err_d1);
    CHARM(err_handler)(err_d1, 1);
    REAL *f = (REAL *)malloc(pnt_d0->local_npoint * sizeof(REAL));
    if (f == NULL)
    {
        fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE);
        exit(CHARM_FAILURE);
    }


    CHARM(shs_point)(pnt_d0, shcs_d1, shcs_d1->nmax, f, err_d1);
    if (err_d1->code == CHARM_SUCCESS)
    {
        fprintf(stderr, "%s",
                "\"charm" CHARM_SUFFIX "_shs_point\" incorrectly "
                "accepted distributed \"charm" CHARM_SUFFIX "_shc\" "
                "and \"charm" CHARM_SUFFIX "_err\" structures and "
                "non-distributed \"charm" CHARM_SUFFIX "_point\" "
                "structure.");
        e += 1;
        CHARM(err_handler)(err_d1, 0);
    }
    CHARM(err_reset)(err_d1);


    CHARM(crd_point_free)(pnt_d0);
    CHARM(shc_free)(shcs_d1);
    free(f);
    }
    /* ..................................................................... */


    /* "charm_shc" is not distributed */
    /* ..................................................................... */
    {
    size_t nlat, nlon;
    CHARM(crd_point_gl_shape)(NMAX_MPI, &nlat, &nlon);
    pnt_d1  = CHARM(mpi_crd_point_gl)(NMAX_MPI, r, nlat * (rank == 0), 0, comm,
                                      err_d1);
    CHARM(err_handler)(err_d1, 1);
    shcs_d0 = CHARM(shc_malloc)(NMAX_MPI, mu, r);
    if (shcs_d0 == NULL)
    {
        fprintf(stderr, "%s", ERR_MSG_SHC);
        exit(CHARM_FAILURE);
    }
    REAL *f = (REAL *)malloc(pnt_d1->local_npoint * sizeof(REAL));
    if (f == NULL)
    {
        fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE);
        exit(CHARM_FAILURE);
    }


    CHARM(shs_point)(pnt_d1, shcs_d0, shcs_d0->nmax, f, err_d1);
    if (err_d1->code == CHARM_SUCCESS)
    {
        fprintf(stderr, "%s",
                "\"charm" CHARM_SUFFIX "_shs_point\" incorrectly "
                "accepted distributed \"charm" CHARM_SUFFIX "_point\" "
                "and \"charm" CHARM_SUFFIX "_err\" structures and "
                "non-distributed \"charm" CHARM_SUFFIX "_shc\" "
                "structure.");
        e += 1;
        CHARM(err_handler)(err_d1, 0);
    }
    CHARM(err_reset)(err_d1);


    CHARM(crd_point_free)(pnt_d1);
    CHARM(shc_free)(shcs_d0);
    free(f);
    }
    /* ..................................................................... */


    /* "charm_err" is not distributed */
    /* ..................................................................... */
    {
    size_t nlat, nlon;
    CHARM(crd_point_gl_shape)(NMAX_MPI, &nlat, &nlon);
    pnt_d1  = CHARM(mpi_crd_point_gl)(NMAX_MPI, r, nlat * (rank == 0), 0, comm,
                                      err_d1);
    CHARM(err_handler)(err_d1, 1);
    unsigned long local_order[2] = {0, NMAX_MPI};
    shcs_d1 = CHARM(mpi_shc_malloc)(NMAX_MPI, mu, r, 1 * (rank == 0),
                                    local_order, comm, err_d1);
    CHARM(err_handler)(err_d1, 1);
    REAL *f = (REAL *)malloc(pnt_d1->local_npoint * sizeof(REAL));
    if (f == NULL)
    {
        fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE);
        exit(CHARM_FAILURE);
    }


    CHARM(shs_point)(pnt_d1, shcs_d1, shcs_d1->nmax, f, err_d0);
    if (err_d0->code == CHARM_SUCCESS)
    {
        fprintf(stderr, "%s",
                "\"charm" CHARM_SUFFIX "_shs_point\" incorrectly "
                "accepted distributed \"charm" CHARM_SUFFIX "_point\" "
                "and \"charm" CHARM_SUFFIX "_shc\" structures and "
                "non-distributed \"charm" CHARM_SUFFIX "_err\" "
                "structure.");
        e += 1;
        CHARM(err_handler)(err_d0, 0);
    }
    CHARM(err_reset)(err_d0);


    CHARM(crd_point_free)(pnt_d1);
    CHARM(shc_free)(shcs_d1);
    free(f);
    }
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(err_free)(err_d1);
    CHARM(err_free)(err_d0);
    /* --------------------------------------------------------------------- */


    return e;
}

