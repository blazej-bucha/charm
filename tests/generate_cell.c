/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/crd/crd_cell_isGrid.h"
#include "../src/crd/crd_cell_isSctr.h"
#include "generate_cell.h"
/* ------------------------------------------------------------------------- */






void CHARM(generate_cell)(CHARM(cell) *grd,
                          REAL r,
                          REAL lat_rng,
                          REAL lon_rng)
{
    const size_t nlat  = grd->nlat;
    const size_t nlon  = grd->nlon;
    const int grd_type = grd->type;


    if (CHARM(crd_cell_isGrid)(grd_type) || CHARM(crd_cell_isSctr)(grd_type))
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
        fprintf(stderr, "Wrong cell type.\n");
        exit(CHARM_FAILURE);
    }


    for (size_t l = 0; l < nlat; l++)
        grd->r[l] = r;


    return;
}
