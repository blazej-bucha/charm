/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdint.h>
#include "../prec.h"
#include "crd_point_isQuadGrid.h"
#include "crd_point_isGLGrid.h"
#include "crd_point_isDH1Grid.h"
#include "crd_point_isDH2Grid.h"
#include "crd_point_quad_equator.h"
/* ------------------------------------------------------------------------- */






/* Returns the index of the equator for the Gauss--Legendre and Driscoll--Healy
 * grids based on the grid type and its maximum harmonic degree.  If if the
 * quadrature grid does not have a zero latitude, returned is "SIZE_MAX". */
size_t CHARM(crd_point_quad_equator)(int grd_type,
                                     unsigned long nmax)
{
    size_t equator;
    if (nmax == CHARM_SHC_NMAX_ERROR)
        equator = SIZE_MAX;
    else if (!CHARM(crd_point_isQuadGrid)(grd_type))
        equator = SIZE_MAX;
    else if (CHARM(crd_point_isGLGrid)(grd_type) && (nmax % 2))
        /* Even number of latitudes in Gauss--Legendre grids means there is no
         * equator. */
        equator = SIZE_MAX;
    else
    {
        size_t nlat = SIZE_MAX;
        size_t nlon = SIZE_MAX;


        if (CHARM(crd_point_isGLGrid)(grd_type))
            CHARM(crd_point_gl_shape)(nmax, &nlat, &nlon);
        else if (CHARM(crd_point_isDH1Grid)(grd_type))
            CHARM(crd_point_dh1_shape)(nmax, &nlat, &nlon);
        else if (CHARM(crd_point_isDH2Grid)(grd_type))
            CHARM(crd_point_dh2_shape)(nmax, &nlat, &nlon);


        equator = nlat / 2;
    }


    return equator;
}
