/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






REAL_SIMD CHARM(misc_polar_optimization_threshold)(unsigned long nmax)
{
    return SET1_R(CHARM_MAX((REAL)CHARM(glob_polar_optimization_a1),
                            CHARM(glob_polar_optimization_a2) * nmax));
}
