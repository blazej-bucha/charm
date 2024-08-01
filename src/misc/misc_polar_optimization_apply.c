/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "misc_polar_optimization_apply.h"
/* ------------------------------------------------------------------------- */






_Bool CHARM(misc_polar_optimization_apply)(unsigned long m,
                                           unsigned long nmax,
                                           REAL_SIMD *sinlat,
                                           size_t nsinlat,
                                           REAL_SIMD threshold)
{
    /* If "glob_polar_optimization_a2" is negative, the polar optimization is
     * not applied */
    if (CHARM(glob_polar_optimization_a2) < PREC(0.0))
        return 0;


    /* Do the test */
    REAL_SIMD nmax_simd = SET1_R(nmax);
    REAL_SIMD m_simd    = SET1_R(m);
    for (size_t i = 0; i < nsinlat; i++)
    {
        MASK2_SIMD test  = LT_R(threshold,
                                SUB_R(m_simd, MUL_R(nmax_simd, sinlat[i])));
        if (!MASK_TRUE_ALL(test))
            return 0;
    }


    return 1;
}
