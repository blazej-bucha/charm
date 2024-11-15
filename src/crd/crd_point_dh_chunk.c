/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../misc/misc_check_radius.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "crd_point_quad_l.h"
#include "crd_point_quad_nlat_north.h"
#include "crd_point_dh_lat_w_chunk.h"
#include "crd_point_dh_chunk.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(crd_point_dh_chunk)(unsigned long nmax,
                                        REAL r,
                                        int dh_type,
                                        size_t local_nlat,
                                        size_t local_0_start,
                                        CHARM(err) *err)
{
    CHARM(point) *dhg = NULL;


    /* Some simple error checks */
    /* --------------------------------------------------------------------- */
    CHARM(misc_check_radius)(r, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto FAILURE;
    }
    /* --------------------------------------------------------------------- */


    /* Initialize a "CHARM(point)" structure based on the "nmax" value. */
    /* --------------------------------------------------------------------- */
    size_t nlat, nlon;
    if (dh_type == CHARM_CRD_POINT_GRID_DH1)
        CHARM(crd_point_dh1_shape)(nmax, &nlat, &nlon);
    else if (dh_type == CHARM_CRD_POINT_GRID_DH2)
        CHARM(crd_point_dh2_shape)(nmax, &nlat, &nlon);
    else
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported value of \"dh_type\".");
        goto FAILURE;
    }


    size_t local_nlat_north;
    if (dh_type == CHARM_CRD_POINT_GRID_DH1)
        local_nlat_north = CHARM(crd_point_quad_nlat_north)(local_nlat,
                                                      local_0_start, nlat,
                                                      CHARM_CRD_POINT_GRID_DH1,
                                                      nmax, err);
    else if (dh_type == CHARM_CRD_POINT_GRID_DH2)
        local_nlat_north = CHARM(crd_point_quad_nlat_north)(local_nlat,
                                                      local_0_start, nlat,
                                                      CHARM_CRD_POINT_GRID_DH2,
                                                      nmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto FAILURE;
    }


    if (dh_type == CHARM_CRD_POINT_GRID_DH1)
    {
        CHARM(crd_point_dh1_shape)(nmax, &nlat, &nlon);
        dhg = CHARM(crd_point_calloc)(CHARM_CRD_POINT_GRID_DH1, local_nlat,
                                      nlon);
    }
    else if (dh_type == CHARM_CRD_POINT_GRID_DH2)
    {
        CHARM(crd_point_dh2_shape)(nmax, &nlat, &nlon);
        dhg = CHARM(crd_point_calloc)(CHARM_CRD_POINT_GRID_DH2, local_nlat,
                                      nlon);
    }
    if (dhg == NULL)
        goto FAILURE;


#if HAVE_MPI
    if (dhg->local_nlat == 0)
        /* This is a valid case, so skip the computation of latitudes,
         * integration weights and spherical radii and continue with
         * longitudes.  */
        goto LONGITUDES;
#endif
    /* --------------------------------------------------------------------- */


    /* Latitudes and integration weights */
    /* --------------------------------------------------------------------- */
    CHARM(crd_point_dh_lat_w_chunk)(dhg, nmax, local_nlat, local_0_start,
                                    local_nlat_north);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto FAILURE;
    }
    /* --------------------------------------------------------------------- */


    /* Spherical radii */
    /* --------------------------------------------------------------------- */
    for (unsigned long i = 0; i < dhg->nlat; i++)
        dhg->r[i] = r;
    /* --------------------------------------------------------------------- */


    /* Longitudes */
    /* --------------------------------------------------------------------- */
#if HAVE_MPI
LONGITUDES:
    ;  /* To avoid declaration after a label */
#endif


    unsigned long L = CHARM(crd_point_quad_l)(nmax);
    REAL c = PREC(0.0);
    if (dh_type == CHARM_CRD_POINT_GRID_DH1)
        c = PI / (REAL)L;
    else if (dh_type == CHARM_CRD_POINT_GRID_DH2)
        c = PI / (REAL)(2 * L);


    for (unsigned long i = 0; i < dhg->nlon; i++)
        dhg->lon[i] = c * (REAL)(i);
    /* --------------------------------------------------------------------- */


EXIT:
    return dhg;


FAILURE:
    CHARM(crd_point_free)(dhg);
    dhg = NULL;


    goto EXIT;
}
