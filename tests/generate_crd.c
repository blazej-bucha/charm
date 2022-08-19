/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(generate_crd)(CHARM(crd) *grd, REAL r, REAL lat_rng, REAL lon_rng)
{
    size_t nlat  = grd->nlat;
    size_t nlon  = grd->nlon;
    int grd_type = grd->type;


    size_t nlat2;
    if ((grd_type == CHARM_CRD_POINTS_GRID) && (nlat > 1))
        /* This ensures that we get symmetric custom point grids */
        nlat2 = nlat - 1;
    else
        nlat2 = nlat;


    if ((grd_type == CHARM_CRD_POINTS_GRID) ||
        (grd_type == CHARM_CRD_POINTS_GRID_GL) ||
        (grd_type == CHARM_CRD_POINTS_GRID_DH1) ||
        (grd_type == CHARM_CRD_POINTS_GRID_DH2) ||
        (grd_type == CHARM_CRD_POINTS_SCATTERED))
    {
        for (size_t l = 0; l < nlat; l++)
            grd->lat[l] = ((REAL)l / (REAL)nlat2) * lat_rng - PI_2;


        for (size_t l = 0; l < nlon; l++)
            grd->lon[l] = ((REAL)l / (REAL)nlon) * lon_rng;
    }
    else if ((grd_type == CHARM_CRD_CELLS_GRID) ||
             (grd_type == CHARM_CRD_CELLS_SCATTERED))
    {
        for (size_t l = 0; l < nlat; l++)
        {
            grd->lat[2 * l]     = ((REAL)l / (REAL)nlat) * lat_rng - PI_2;
            grd->lat[2 * l + 1] = ((REAL)(l + 1) / (REAL)nlat) * lat_rng -
                                  PI_2;
        }


        for (size_t l = 0; l < nlon; l++)
        {
            grd->lon[2 * l]     = ((REAL)l / (REAL)nlon) * lon_rng;
            grd->lon[2 * l + 1] = ((REAL)(l + 1) / (REAL)nlon) * lon_rng;
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
