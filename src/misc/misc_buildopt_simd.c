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
#elif HAVE_NEON
    return BUILDOPT_SIMD_NEON;
#elif HAVE_SSE41
    return BUILDOPT_SIMD_SSE41;
#else
    return BUILDOPT_SIMD_NONE;
#endif
}
