/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "crd_point_dh_lats_weights.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(crd_point_dh2)(unsigned long nmax, REAL r)
{
    /* Initialize a "CHARM(point)" structure based on the "nmax" value. */
    /* --------------------------------------------------------------------- */
    unsigned long L = nmax + 1;


    CHARM(point) *dhg = CHARM(crd_point_calloc)(CHARM_CRD_POINT_GRID_DH2,
                                                2 * L, 4 * L);
    if (dhg == NULL)
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Latitudes and integration weights */
    /* --------------------------------------------------------------------- */
    CHARM(crd_point_dh_lats_weights)(dhg, nmax);
    /* --------------------------------------------------------------------- */


    /* Longitudes */
    /* --------------------------------------------------------------------- */
    REAL c = PI / (REAL)(2 * L);


#if CHARM_OPENMP
#pragma omp parallel for default(none) shared(dhg, L, c)
#endif
    for (unsigned long i = 0; i < (4 * L); i++)
        dhg->lon[i] = c * (REAL)(i);
    /* --------------------------------------------------------------------- */


    /* Spherical radii */
    /* --------------------------------------------------------------------- */
#if CHARM_OPENMP
#pragma omp parallel for default(none) shared(dhg, L, r)
#endif
    for (unsigned long i = 0; i < (2 * L); i++)
        dhg->r[i] = r;
    /* --------------------------------------------------------------------- */


    return dhg;
}
