/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "misc_buildopt.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_buildopt_simd)(void)
{
#if HAVE_AVX_INSTRUCTIONS
    return BUILDOPT_SIMD_AVX;
#elif HAVE_AVX2_INSTRUCTIONS
    return BUILDOPT_SIMD_AVX2;
#elif HAVE_AVX512F_INSTRUCTIONS
    return BUILDOPT_SIMD_AVX512;
#endif
#if !(HAVE_AVX_INSTRUCTIONS) && !(HAVE_AVX2_INSTRUCTIONS) && \
    !(HAVE_AVX512F_INSTRUCTIONS)
    return BUILDOPT_SIMD_NONE;
#endif
}
