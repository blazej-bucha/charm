/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/crd/crd_point_isGrid.h"
#include "../src/crd/crd_point_isSctr.h"
#include "generate_point.h"
/* ------------------------------------------------------------------------- */






void CHARM(generate_point)(CHARM(point) *pnt,
                           REAL r,
                           REAL lat_rng,
                           REAL lon_rng)
{
#if HAVE_MPI
    const size_t nlat  = pnt->local_nlat;
    const size_t nlon  = pnt->local_nlon;
#else
    const size_t nlat  = pnt->nlat;
    const size_t nlon  = pnt->nlon;
#endif
    const int grd_type = pnt->type;


    size_t nlat2;
    if (CHARM(crd_point_isGrid)(grd_type) && (nlat > 1))
        /* This ensures that we get symmetric custom point grids */
        nlat2 = nlat - 1;
    else
        nlat2 = nlat;


    if (CHARM(crd_point_isGrid)(grd_type) || CHARM(crd_point_isSctr)(grd_type))
    {
        for (size_t l = 0; l < nlat; l++)
            pnt->lat[l] = PI_2 - ((REAL)l / (REAL)nlat2) * lat_rng;


        for (size_t l = 0; l < nlon; l++)
            pnt->lon[l] = ((REAL)l / (REAL)nlon) * lon_rng;
    }
    else
    {
        fprintf(stderr, "Wrong point type.\n");
        exit(CHARM_FAILURE);
    }


    for (size_t l = 0; l < nlat; l++)
        pnt->r[l] = r;


    return;
}
