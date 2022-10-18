/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(generate_cell)(CHARM(cell) *grd, REAL r, REAL lat_rng, REAL lon_rng)
{
    size_t nlat  = grd->nlat;
    size_t nlon  = grd->nlon;
    int grd_type = grd->type;


    if ((grd_type == CHARM_CRD_CELLS_GRID) ||
        (grd_type == CHARM_CRD_CELLS_SCATTERED))
    {
        for (size_t l = 0; l < nlat; l++)
        {
            grd->latmin[l] = PI_2 - ((REAL)(l + 1) / (REAL)nlat) * lat_rng;
            grd->latmax[l] = PI_2 - ((REAL)l / (REAL)nlat) * lat_rng;
        }


        for (size_t l = 0; l < nlon; l++)
        {
            grd->lonmin[l] = ((REAL)l / (REAL)nlon) * lon_rng;
            grd->lonmax[l] = ((REAL)(l + 1) / (REAL)nlon) * lon_rng;
        }
    }
    else
    {
        fprintf(stderr, "Wrong grid type.\n");
        exit(1);
    }


    for (size_t l = 0; l < nlat; l++)
        grd->r[l] = r;


    return;
}
