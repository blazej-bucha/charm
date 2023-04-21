/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(generate_point)(CHARM(point) *grd, REAL r, REAL lat_rng,
                           REAL lon_rng)
{
    size_t nlat  = grd->nlat;
    size_t nlon  = grd->nlon;
    int grd_type = grd->type;


    size_t nlat2;
    if ((grd_type == CHARM_CRD_POINT_GRID) && (nlat > 1))
        /* This ensures that we get symmetric custom point grids */
        nlat2 = nlat - 1;
    else
        nlat2 = nlat;


    if ((grd_type == CHARM_CRD_POINT_GRID) ||
        (grd_type == CHARM_CRD_POINT_SCATTERED))
    {
        for (size_t l = 0; l < nlat; l++)
            grd->lat[l] = PI_2 - ((REAL)l / (REAL)nlat2) * lat_rng;


        for (size_t l = 0; l < nlon; l++)
            grd->lon[l] = ((REAL)l / (REAL)nlon) * lon_rng;
    }
    else
    {
        fprintf(stderr, "Wrong grid type.\n");
        exit(CHARM_FAILURE);
    }


    for (size_t l = 0; l < nlat; l++)
        grd->r[l] = r;


    return;
}
