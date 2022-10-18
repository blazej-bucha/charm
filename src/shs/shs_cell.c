/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shs_cell_grd.h"
#include "shs_cell_sctr.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell)(const CHARM(cell) *cell, const CHARM(shc) *shcs,
                     unsigned long nmax, REAL *f, CHARM(err) *err)
{
    /* Some trivial initial error checks */
    /* --------------------------------------------------------------------- */
    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Maximum harmonic degree of the synthesis "
                       "(\"nmax\") cannot be larger than "
                       "maximum harmonic degree of spherical harmonic "
                       "coefficients (\"shcs->nmax\").");
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Now do the synthesis */
    /* --------------------------------------------------------------------- */
    if (cell->type == CHARM_CRD_CELLS_SCATTERED)
    {
        if (cell->nlat != cell->nlon)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The number of latitudes and longitudes in the "
                           "\"cell\" structure must be the same in order to "
                           "perform cell-wise spherical harmonic synthesis.");
            return;
        }


        /* Point-wise synthesis */
        CHARM(shs_cell_sctr)(cell, shcs, nmax, f, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
    else if (cell->type == CHARM_CRD_CELLS_GRID)
    {
        /* Grid-wise synthesis */
        CHARM(shs_cell_grd)(cell, shcs, nmax, f, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
    else
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported \"cell->type\" for spherical harmonic "
                       "synthesis of block-mean values in cells.");
        return;
    }
    /* --------------------------------------------------------------------- */






    return;
}
