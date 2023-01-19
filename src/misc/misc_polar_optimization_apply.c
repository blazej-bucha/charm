/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






_Bool CHARM(misc_polar_optimization_apply)(unsigned long m,
                                           unsigned long nmax,
                                           REAL_SIMD sinlat,
                                           REAL_SIMD threshold)
{
    /* If "glob_polar_optimization_a2" is negative, the polar optimization is
     * not applied */
    if (CHARM(glob_polar_optimization_a2) < PREC(0.0))
        return 0;


    /* Do the test */
    REAL_SIMD nmax_simd = SET1_R(nmax);
    REAL_SIMD m_simd    = SET1_R(m);
    MASK2_SIMD test     = LT_R(threshold,
                               SUB_R(m_simd, MUL_R(nmax_simd, sinlat)));


    /* Apply the polar optimization only if "test" is true for all "sinlat". */
    if (MOVEMASK((test)) == SIMD_TRUE)
        return 1;
    else
        return 0;
}
