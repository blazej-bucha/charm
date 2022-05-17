/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "crd_dh_lats_weights.h"
/* ------------------------------------------------------------------------- */







CHARM(crd) *CHARM(crd_dh1)(unsigned long nmax, REAL r)
{
    /* Initialize a "CHARM(crd)" structure based on the "nmax" value. */
    /* --------------------------------------------------------------------- */
    unsigned long L = nmax + 1;


    CHARM(crd) *dhg = CHARM(crd_init)(CHARM_CRD_POINTS_GRID_DH1, 2 * L, 2 * L);
    if (dhg == NULL)
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Latitudes and integration weights */
    /* --------------------------------------------------------------------- */
    CHARM(crd_dh_lats_weights)(dhg, nmax);
    /* --------------------------------------------------------------------- */


    /* Longitudes */
    /* --------------------------------------------------------------------- */
    REAL c = PI / (REAL)L;


#if CHARM_PARALLEL
#pragma omp parallel for default(none) shared(dhg, L, c)
#endif
    for (unsigned long i = 0; i < (2 * L); i++)
        dhg->lon[i] = c * (REAL)(i);
    /* --------------------------------------------------------------------- */


    /* Spherical radii */
    /* --------------------------------------------------------------------- */
#if CHARM_PARALLEL
#pragma omp parallel for default(none) shared(dhg, L, r)
#endif
    for (unsigned long i = 0; i < (2 * L); i++)
        dhg->r[i] = r;
    /* --------------------------------------------------------------------- */


    return dhg;
}
