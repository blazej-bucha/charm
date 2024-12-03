/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_isGLGrid.h"
#include "crd_point_isDHGrid.h"
#include "crd_point_quad_get_nmax_from_nlat.h"
/* ------------------------------------------------------------------------- */






/* Returns the maximum harmonic degree of quadrature grids based on the grid
 * type and its number of latitudes.  For unsupported point types,
 * "CHARM_SHC_NMAX_ERROR" is returned. */
size_t CHARM(crd_point_quad_get_nmax_from_nlat)(int grd_type,
                                                size_t nlat)
{
    size_t nmax;
    if (CHARM(crd_point_isGLGrid)(grd_type))
        nmax = nlat - 1;
    else if (CHARM(crd_point_isDHGrid)(grd_type))
        nmax = (nlat - 2) / 2;
    else
        nmax = CHARM_SHC_NMAX_ERROR;


    return nmax;
}
