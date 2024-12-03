/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isQuadGrid.h"
#include "crd_point_isGLGrid.h"
#include "crd_point_isDHGrid.h"
#include "crd_point_isGrid.h"
#include "crd_cell_isGrid.h"
#include "crd_grd_check_symm.h"
/* ------------------------------------------------------------------------- */






/* An internal function to determine whether or not to employ the symmetry
 * property of Legendre functions. */
void CHARM(crd_grd_check_symm)(size_t ipv,
                               size_t v,
                               size_t local_0_start,
                               size_t equator,
                               int grd_type,
                               size_t nlatdo,
                               _Bool symm,
                               _Bool even,
                               REAL *symmv,
                               REAL *latsinv)
{
    if (ipv >= nlatdo)
    {
        /* We are outside the bounds of the input latitudes. */
        latsinv[v] = 0;
        symmv[v]   = 0;
        return;
    }
    else
        latsinv[v] = 1;


    if (CHARM(crd_point_isQuadGrid)(grd_type))
    {
        /* Get the global index of the latitude */
        size_t lpipv = local_0_start + ipv;


        if (CHARM(crd_point_isDHGrid)(grd_type) && (lpipv == 0))
            /* The first latitude of Driscoll--Healy grids (the north pole)
             * does not have its negative counterpart (the south pole), so no
             * symmetry. */
            symmv[v] = 0;
        else if (lpipv < equator)
            symmv[v] = 1;
        else /* if (lpipv == equator) */
            symmv[v] = 0;
    }
    else if (CHARM(crd_point_isGrid)(grd_type) ||
             CHARM(crd_cell_isGrid)(grd_type))
    {
        if (!symm)
            symmv[v] = 0;
        else if (symm && !even && (ipv == (nlatdo - 1)))
            symmv[v] = 0;
        else
            symmv[v] = 1;
    }


    return;
}

