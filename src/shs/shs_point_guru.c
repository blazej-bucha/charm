/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shs_point_grd.h"
#include "shs_point_sctr.h"
#include "shs_check_single_derivative.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../crd/crd_point_isSctr.h"
#include "../crd/crd_point_isGrid.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_point_guru)(const CHARM(point) *pnt,
                           const CHARM(shc) *shcs,
                           unsigned long nmax,
                           unsigned dr,
                           unsigned dlat,
                           unsigned dlon,
                           REAL *f,
                           CHARM(err) *err)
{
    /* Some trivial initial error checks */
    /* --------------------------------------------------------------------- */
    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Maximum harmonic degree of the synthesis (\"nmax\") "
                       "cannot be larger than maximum harmonic degree of "
                       "spherical harmonic coefficients (\"shcs->nmax\").");
        return;
    }


    /* Importantly, "dr", "dlat" and "dlon" will be cast to "int" when passing
     * them to "shs_point_*" functions below, which may take also negative ints
     * as special values.  However, this function assumes non-negative ints,
     * hence the "unsigned int" data type.  Now check for correct values of
     * "dr", "dlat" and "dlon". */
    CHARM(shs_check_single_derivative)(dr, dlat, dlon, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Do nothing if the total number of points in "pnt" is zero, which is
     * a valid case */
    /* --------------------------------------------------------------------- */
    if (pnt->npoint == 0)
        return;
    /* --------------------------------------------------------------------- */






    /* Now do the synthesis */
    /* --------------------------------------------------------------------- */
    if (CHARM(crd_point_isSctr)(pnt->type))
    {
        if (pnt->nlat != pnt->nlon)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The number of latitudes and longitudes in the "
                           "\"pnt\" structure must be the same to "
                           "perform point-wise spherical harmonic synthesis.");
            return;
        }


        /* Point-wise synthesis */
        CHARM(shs_point_sctr)(pnt, shcs, nmax, dr, dlat, dlon, &f, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
    else if (CHARM(crd_point_isGrid)(pnt->type))
    {
        /* Grid-wise synthesis */
        CHARM(shs_point_grd)(pnt, shcs, nmax, dr, dlat, dlon, &f, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
    else
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported \"pnt->type\" for spherical harmonic "
                       "synthesis of point values.");
        return;
    }
    /* --------------------------------------------------------------------- */






    return;
}
