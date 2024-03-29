/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shs_point_grd.h"
#include "shs_point_sctr.h"
#include "shs_point_gradn.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../crd/crd_point_isSctr.h"
#include "../crd/crd_point_isGrid.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_point)(const CHARM(point) *pnt,
                      const CHARM(shc) *shcs,
                      unsigned long nmax,
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
        CHARM(shs_point_sctr)(pnt, shcs, nmax, GRAD_0, GRAD_0, GRAD_0, &f,
                              err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
    else if (CHARM(crd_point_isGrid)(pnt->type))
    {
        /* Grid-wise synthesis */
        CHARM(shs_point_grd)(pnt, shcs, nmax, GRAD_0, GRAD_0, GRAD_0, &f, err);
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
