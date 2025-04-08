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
#include "../src/shc/shc_reset_coeffs.h"
#include "../src/err/err_is_null_ptr.h"
#include "../src/crd/crd_point_isGLGrid.h"
#include "../src/crd/crd_point_isDH1Grid.h"
#include "../src/crd/crd_point_isDH2Grid.h"
#include "parameters.h"
#include "error_messages.h"
#include "mpi_shc_distribute.h"
#include "mpi_crd_point_distribute.h"
#include "cmp_arrays.h"
#include "check_mpi_sha_point.h"
/* ------------------------------------------------------------------------- */






/* ------------------------------------------------------------------------- */
#undef NPOINT_TYPE
#define NPOINT_TYPE (3)
/* ------------------------------------------------------------------------- */






long int check_mpi_sha_point(void)
{
    /* --------------------------------------------------------------------- */
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    long int e = 0;


    const REAL mu = PREC(1.1);
    const REAL r  = PREC(1.2);
    const REAL dr = PREC(0.1);


    CHARM(point) *grd = NULL;


    CHARM(shc) *shcs_d0     = NULL;
    CHARM(shc) *shcs_d1     = NULL;
    CHARM(shc) *shcs_d1_ret = NULL;


    int point_types[NPOINT_TYPE] = {CHARM_CRD_POINT_GRID_GL,
                                    CHARM_CRD_POINT_GRID_DH1,
                                    CHARM_CRD_POINT_GRID_DH2};
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


    /* --------------------------------------------------------------------- */
    unsigned block_a = CHARM(glob_sha_block_lat_multiplier);


    for (int t = 0; t < NPOINT_TYPE; t++)
    {
    int type = point_types[t];


    for (unsigned b = 0; b < 3; b++)
    {
    CHARM(glob_sha_block_lat_multiplier) = block_a + b;


    for (int x = 0; x < 2; x++)
    {
    /* If "x == 0", all spherical harmonic coefficients will be distributed on
     * one MPI process only.  If "x == 1", the coefficients are distributed
     * across all MPI processes. */
    for (unsigned long nmax = NMAX_MPI; nmax <= NMAX_MPI + 3; nmax++)
    {


    size_t nlat, nlon, local_nlat, local_0_start;
    if (CHARM(crd_point_isGLGrid)(type))
    {
        CHARM(crd_point_gl_shape)(nmax, &nlat, &nlon);
        mpi_crd_point_distribute(nlat, type, rank, size, &local_nlat,
                                 &local_0_start);
        grd = CHARM(mpi_crd_point_gl)(nmax, r + dr, local_nlat, local_0_start,
                                      comm, err_d1);
        CHARM(err_handler)(err_d1, 1);
    }
    else if (CHARM(crd_point_isDH1Grid)(type))
    {
        CHARM(crd_point_dh1_shape)(nmax, &nlat, &nlon);
        mpi_crd_point_distribute(nlat, type, rank, size, &local_nlat,
                                 &local_0_start);
        grd = CHARM(mpi_crd_point_dh1)(nmax, r + dr, local_nlat, local_0_start,
                                      comm, err_d1);
        CHARM(err_handler)(err_d1, 1);
    }
    else if (CHARM(crd_point_isDH2Grid)(type))
    {
        CHARM(crd_point_dh2_shape)(nmax, &nlat, &nlon);
        mpi_crd_point_distribute(nlat, type, rank, size, &local_nlat,
                                 &local_0_start);
        grd = CHARM(mpi_crd_point_dh2)(nmax, r + dr, local_nlat, local_0_start,
                                      comm, err_d1);
        CHARM(err_handler)(err_d1, 1);
    }


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


    /* Create some fake spherical harmonic coefficients.  See
     * "check_mpi_shs_point" for details. */
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


    shcs_d1_ret = CHARM(mpi_shc_calloc)(nmax, mu, r, nchunk_tmp,
                                        local_order, comm, err_d1);
    CHARM(err_handler)(err_d1, 1);


    free(local_order);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    REAL *f = (REAL *)malloc(grd->local_npoint * sizeof(REAL));
    if (f == NULL)
    {
        fprintf(stderr, "%s", CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    CHARM(shs_point)(grd, shcs_d1, shcs_d1->nmax, f, err_d1);
    CHARM(err_handler)(err_d1, 1);


    CHARM(sha_point)(grd, f, shcs_d1_ret->nmax, shcs_d1_ret, err_d1);
    CHARM(err_handler)(err_d1, 1);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    /* Before calling "cmp_arrays", keep in mind that "shcs_d1_ret->local_nc ==
     * shcs_d1->local_nc", etc..  Just to be sure, we check this. */
    if ((shcs_d1_ret->local_nc != shcs_d1->local_nc) ||
        (shcs_d1_ret->local_ns != shcs_d1->local_ns))
    {
        fprintf(stderr, "%s",
                "The number of local spherical harmonic coefficients "
                "of a distributed structure does not match the number "
                "of coefficients of a non-distributed structure.\n");
        exit(CHARM_FAILURE);
    }


    if (shcs_d1_ret->local_nchunk > 0)
    {
        e += cmp_arrays(shcs_d1_ret->c[shcs_d1_ret->local_order[0]],
                        shcs_d1->c[shcs_d1->local_order[0]],
                        shcs_d1_ret->local_nc,
                        CHARM(glob_threshold2));
        e += cmp_arrays(shcs_d1_ret->s[shcs_d1_ret->local_order[0]],
                        shcs_d1->s[shcs_d1->local_order[0]],
                        shcs_d1_ret->local_ns,
                        CHARM(glob_threshold2));
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    free(f);
    CHARM(shc_free)(shcs_d0);
    CHARM(shc_free)(shcs_d1);
    CHARM(shc_free)(shcs_d1_ret);
    /* --------------------------------------------------------------------- */
    }


    CHARM(crd_point_free)(grd);
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
    CHARM(point) *pnt_d0  = CHARM(crd_point_gl)(NMAX_MPI, r);
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


    CHARM(sha_point)(pnt_d0, f, shcs_d1->nmax, shcs_d1, err_d1);
    if (err_d1->code == CHARM_SUCCESS)
    {
        fprintf(stderr, "%s",
                "\"charm" CHARM_SUFFIX "_sha_point\" incorrectly "
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
    CHARM(point) *pnt_d1  = CHARM(mpi_crd_point_gl)(NMAX_MPI, r,
                                                    nlat * (rank == 0),
                                                    0, comm, err_d1);
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


    CHARM(sha_point)(pnt_d1, f, shcs_d0->nmax, shcs_d0, err_d1);
    if (err_d1->code == CHARM_SUCCESS)
    {
        fprintf(stderr, "%s",
                "\"charm" CHARM_SUFFIX "_sha_point\" incorrectly "
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
    CHARM(point) *pnt_d1  = CHARM(mpi_crd_point_gl)(NMAX_MPI, r,
                                                    nlat * (rank == 0),
                                                    0, comm, err_d1);
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


    CHARM(sha_point)(pnt_d1, f, shcs_d1->nmax, shcs_d1, err_d0);
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

