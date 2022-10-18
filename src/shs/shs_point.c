/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shs_point_grd.h"
#include "shs_point_sctr.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_point)(const CHARM(point) *pnt, const CHARM(shc) *shcs,
                      unsigned long nmax, REAL *f, CHARM(err) *err)
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
    /* --------------------------------------------------------------------- */






    /* Now do the synthesis */
    /* --------------------------------------------------------------------- */
    if (pnt->type == CHARM_CRD_POINTS_SCATTERED)
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
        CHARM(shs_point_sctr)(pnt, shcs, nmax, f, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
    else if ((pnt->type == CHARM_CRD_POINTS_GRID) ||
             (pnt->type == CHARM_CRD_POINTS_GRID_GL) ||
             (pnt->type == CHARM_CRD_POINTS_GRID_DH1) ||
             (pnt->type == CHARM_CRD_POINTS_GRID_DH2))
    {
        /* Grid-wise synthesis */
        CHARM(shs_point_grd)(pnt, shcs, nmax, f, err);
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
