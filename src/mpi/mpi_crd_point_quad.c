/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if HAVE_MPI
#   include <mpi.h>
#else
#   error "mpi.h required"
#endif
#include "../prec.h"
#include "../crd/crd_point_isGLGrid.h"
#include "../crd/crd_point_isDHGrid.h"
#include "../crd/crd_point_gl_chunk.h"
#include "../crd/crd_point_dh1_chunk.h"
#include "../crd/crd_point_dh2_chunk.h"
#include "../crd/crd_point_issymm.h"
#include "../misc/misc_is_nearly_equal.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "mpi_crd_point_check_struct.h"
#include "mpi_crd_point_local2distributed.h"
#include "mpi_size_t.h"
#include "mpi_allequal.h"
#include "mpi_crd_point_quad.h"
#include "mpi_err_isdistributed.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(mpi_crd_point_quad)(unsigned long nmax,
                                        REAL r,
                                        size_t local_nlat,
                                        size_t local_0_start,
                                        MPI_Comm comm,
                                        void (*quad_shape)(unsigned long,
                                                           size_t *,
                                                           size_t *),
                                     CHARM(point) *(*quad_chunk)(unsigned long,
                                                                 REAL,
                                                                 size_t,
                                                                 size_t,
                                                                 CHARM(err) *),
                                        CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    char err_msg[CHARM_ERR_MAX_MSG];
    CHARM(point) *pnt = NULL;
    int grd_type = 0;


    if (!CHARM(mpi_err_isdistributed)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto BARRIER;
    }
    /* --------------------------------------------------------------------- */


    /* Identify grid type */
    /* --------------------------------------------------------------------- */
    if ((quad_shape == CHARM(crd_point_gl_shape)) &&
        (quad_chunk == CHARM(crd_point_gl_chunk)))
        grd_type = CHARM_CRD_POINT_GRID_GL;
    else if ((quad_shape == CHARM(crd_point_dh1_shape)) &&
             (quad_chunk == CHARM(crd_point_dh1_chunk)))
        grd_type = CHARM_CRD_POINT_GRID_DH1;
    else if ((quad_shape == CHARM(crd_point_dh2_shape)) &&
             (quad_chunk == CHARM(crd_point_dh2_chunk)))
        grd_type = CHARM_CRD_POINT_GRID_DH2;
    else
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong value of \"quad_shape\" and/or \"quad_chunk\".");


BARRIER:
    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* "nmax" and "r" across all processes must be equal */
    /* --------------------------------------------------------------------- */
    /* ..................................................................... */
    _Bool ret = CHARM(mpi_allequal_ulong)(nmax, comm, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM(mpi_err_isempty)(err))
            goto FAILURE;


    if (!ret)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Maximum harmonic degree \"nmax\" must be equal "
                       "for all processes in \"comm\".");


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* ..................................................................... */


    /* ..................................................................... */
    ret = CHARM(mpi_allequal_real)(r, comm, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;


    if (!ret)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Spherical radius \"r\" must be equal for all "
                       "processes in \"comm\".");


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */


    /* Check that the sum of "local_nlat" across all processes matches the
     * number of latitudes in the grid for "nmax" */
    /* --------------------------------------------------------------------- */
    size_t nlat, nlon;
    quad_shape(nmax, &nlat, &nlon);


    size_t local_nlat_sum;
    MPI_Allreduce(&local_nlat, &local_nlat_sum, 1, CHARM_MPI_SIZE_T, MPI_SUM,
                  comm);
    if (local_nlat_sum != nlat)
    {
        sprintf(err_msg, "The sum of \"local_nlat\" is \"%zu\" which does not "
                         "match the number of latitudes \"%zu\" of this "
                         "grid for \"nmax = %lu\".",
                         local_nlat_sum, nlat, nmax);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
    }


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* Check "local_nlat" */
    /* --------------------------------------------------------------------- */
    if (CHARM(crd_point_isGLGrid)(grd_type))
    {
        if (nmax % 2)  /* Odd "nmax" */
        {
            /* For odd "nmax", all "local_nlat" must be even */
            _Bool local_even = !(local_nlat % 2);
            _Bool even;
            MPI_Allreduce(&local_even, &even, 1, MPI_C_BOOL, MPI_LAND, comm);


            if (!even)
            {
                sprintf(err_msg, "For Gauss--Legendre grids and odd "
                                 "\"nmax = %lu\", all processes in \"comm\" "
                                 "must have even \"local_nlat\".", nmax);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG, err_msg);
            }
        }
        else  /* Even "nmax" */
        {
            /* For even "nmax", "local_nlat" must be odd on one process exactly
             * */

            /* "int" is used intentionally instead of "_Bool" to allow simple
             * sum across all processes */
            int local_odd = local_nlat % 2;
            int odd;
            MPI_Allreduce(&local_odd, &odd, 1, MPI_INT, MPI_SUM, comm);


            if (odd != 1)
            {
                sprintf(err_msg, "For Gauss--Legendre grids and even "
                                 "\"nmax = %lu\", all but one process in "
                                 "\"comm\" must have even \"local_nlat\".",
                                 nmax);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG, err_msg);
            }
        }
    }
    else if (CHARM(crd_point_isDHGrid)(grd_type))
    {
        /* Check if one process holds entire DH grid */
        size_t local_nlat_sum;
        MPI_Allreduce(&local_nlat, &local_nlat_sum, 1, CHARM_MPI_SIZE_T,
                      MPI_SUM, comm);


        int local_nlat_nonzero = local_nlat != 0;
        int local_nlat_nonzero_sum;
        MPI_Allreduce(&local_nlat_nonzero, &local_nlat_nonzero_sum, 1,
                      MPI_INT, MPI_SUM, comm);


        int whole_grd = (local_nlat_sum == nlat) && (local_nlat_nonzero == 1);
        int whole_grd_sum;
        MPI_Allreduce(&whole_grd, &whole_grd_sum, 1, MPI_INT, MPI_SUM, comm);


        if (whole_grd_sum)
        {
            /* One MPI process holds entire DH grid */
            if (local_nlat == nlat)
            {
                if (local_0_start != 0)
                {
                    sprintf(err_msg, "It seems that one MPI process is asked "
                                     "to hold the entire Driscoll--Healy "
                                     "grid.  On one process, "
                                     "it must then hold \"local_0_start = 0\" "
                                     "and \"local_nlat = %zu\", the latter of "
                                     "which is the total number of latitudes "
                                     "for \"nmax = %zu\".", nlat, nmax);
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EFUNCARG, err_msg);
                }
            }
        }
        else
        {
            /* DH grid is stored across MPI processes */


            /* Two processes must have odd "local_nlat", one that stores the
             * north pole and the other that stores the equator. */
            int local_nlat_odd = local_nlat % 2;
            int local_nlat_odd_sum;
            MPI_Allreduce(&local_nlat_odd, &local_nlat_odd_sum, 1,
                          MPI_INT, MPI_SUM, comm);
            if (local_nlat_odd_sum != 2)
            {
                sprintf(err_msg, "For Driscoll--Healy grids distributed "
                                 "across more than one MPI process, "
                                 "\"local_nlat\" must be odd on \"2\" "
                                 "processes and not on \"%d\" "
                                 "processes.", local_nlat_odd_sum);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG, err_msg);
            }
        }
    }


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* Sanity check for "local_0_start" */
    /* --------------------------------------------------------------------- */
    /* Index of the last non-negative latitude */
    size_t last_nonneg_lat = 0;
    if (CHARM(crd_point_isGLGrid)(grd_type))
        last_nonneg_lat = (nlat - 1) / 2;
    else if (CHARM(crd_point_isDHGrid)(grd_type))
        last_nonneg_lat = nlat / 2;



    if ((local_nlat > 0) && (local_0_start > last_nonneg_lat))
    {
        sprintf(err_msg, "For quadrature grids, \"local_0_start\" "
                         "cannot point to a negative latitude.  "
                         "In this case, \"local_0_start\" is \"%zu\" but must "
                         "not exceed \"%zu\", which is the index of "
                         "the last non-negative latitude for this grid type "
                         "and \"nmax = %lu\".",
                         local_0_start, last_nonneg_lat, nmax);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
    }


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    pnt = quad_chunk(nmax, r, local_nlat, local_0_start, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    /* Barrier before a collective call */
    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;


    CHARM(mpi_crd_point_local2distributed)(pnt, local_0_start, comm);
    /* --------------------------------------------------------------------- */


    /* Check the distribution of the latitudinal chunks */
    /* ----------------------------------------------------------------- */
    CHARM(mpi_crd_point_check_struct)(pnt, 1, err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;


    /* The grid must be symmetric with respect to the equator, so let's check
     * this.  We do not call here "mpi_crd_point_issymm", as it intentionally
     * returns "1" if "local_nlat == 0" at some process (provided that "pnt" is
     * symmetric for other processes).  This case must be treated separatelly,
     * hence "crd_point_check_issymm". */
    _Bool symm = CHARM(crd_point_issymm)(pnt);
    if (!symm && (pnt->local_nlat > 1))
    {
        if (CHARM(crd_point_isGLGrid)(pnt->type))
        {
            /* "pnt" must be symmetric if there is more than one latitude. */
            sprintf(err_msg, "The chunk with \"local_0_start = %zu\" and "
                             "\"local_nlat = %zu\" lacks equatorial symmetry.",
                             local_0_start, local_nlat);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
        }
        else if (CHARM(crd_point_isDHGrid)(pnt->type))
        {
            if (pnt->local_0_start != 0)
            {
                /* All latitudinal chunks except the first one with
                 * "pnt->local_0_start" and "pnt->local_nlat > 1" must be
                 * symmetric. */
                sprintf(err_msg, "The chunk with \"local_0_start = %zu\" and "
                                 "\"local_nlat = %zu\" lacks equatorial "
                                 "symmetry.",
                                 local_0_start, local_nlat);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFUNCARG, err_msg);
            }
        }
    }


    if (!CHARM(mpi_err_isempty)(err))
        goto FAILURE;
    /* ----------------------------------------------------------------- */


EXIT:
    return pnt;


FAILURE:
    CHARM(crd_point_free)(pnt);
    pnt = NULL;
    goto EXIT;
}
