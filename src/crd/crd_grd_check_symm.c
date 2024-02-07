/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_grd_check_symm.h"
/* ------------------------------------------------------------------------- */






/* An internal function to determine whether or not the symmetry property of
 * Legendre functions. */
void CHARM(crd_grd_check_symm)(size_t ipv, size_t v, int grd_type,
                               size_t nlatdo, _Bool symm, _Bool even,
                               REAL *symmv, REAL *latsinv)
{
    if ((grd_type == CHARM_CRD_POINT_GRID_DH1) ||
        (grd_type == CHARM_CRD_POINT_GRID_DH2))
    {
        if (ipv == 0)
        {
            /* Do not apply the symmetry property at the north pole
             * of the Driscoll--Healy grids, as they do not have
             * their negative counterpart (the south pole). */
            symmv[v]   = 0;
            latsinv[v] = 1;
        }
        else if (ipv < (nlatdo - 1))
        {
            /* Any latitude between the north pole and the equator
             * of the DH grids.  */
            symmv[v]   = 1;
            latsinv[v] = 1;
        }
        else if (ipv == (nlatdo - 1))
        {
            /* Equator of the DH grids, so no symmetry */
            symmv[v]   = 0;
            latsinv[v] = 1;
        }
        else
        {
            symmv[v]   = 0;
            latsinv[v] = 0;
        }
    }
    else if (grd_type == CHARM_CRD_POINT_GRID_GL)
    {
        if (ipv < (nlatdo + even - 1))
        {
            /* Any latitude up to the equator (excluding).  */
            symmv[v]   = 1;
            latsinv[v] = 1;
        }
        else if (ipv == (nlatdo - 1))
        {
            /* Equator of the GL grid, so no symmetry */
            symmv[v]   = 0;
            latsinv[v] = 1;
        }
        else
        {
            symmv[v]   = 0;
            latsinv[v] = 0;
        }
    }
    else if ((grd_type == CHARM_CRD_POINT_GRID) ||
             (grd_type == CHARM_CRD_CELL_GRID))
    {
        if (symm == 0)
            symmv[v] = 0;
        else if (symm == 1 && even == 0 && (ipv == (nlatdo - 1)))
            symmv[v] = 0;
        else
            symmv[v] = 1;


        if (ipv < nlatdo)
            latsinv[v] = 1;
        else
            latsinv[v] = 0;
    }
    else
    {
        symmv[v] = latsinv[v] = 0;
    }


    return;
}

