/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_quad_l.h"
/* ------------------------------------------------------------------------- */







/* Returns the "L" parameter determining the number of latitudes and longitudes
 * of the GL and DH grids based on "nmax" */
unsigned long CHARM(crd_point_quad_l)(unsigned long nmax)
{
    return nmax + 1;
}
