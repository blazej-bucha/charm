/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_buildopt_simd)(void)
{
#if HAVE_AVX
    return BUILDOPT_SIMD_AVX;
#elif HAVE_AVX2
    return BUILDOPT_SIMD_AVX2;
#elif HAVE_AVX512F
    return BUILDOPT_SIMD_AVX512;
#endif
#if !(HAVE_AVX) && !(HAVE_AVX2) && !(HAVE_AVX512F)
    return BUILDOPT_SIMD_NONE;
#endif
}
