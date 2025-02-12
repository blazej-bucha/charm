/* This header file is not a part of API. */


#ifndef __MISC_BUILDOPT_H__
#define __MISC_BUILDOPT_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/* Signalizes the feature is available.  Must be non-zero. */
#undef BUILDOPT_OMP_CHARM
#define BUILDOPT_OMP_CHARM 1


/* Signalizes the feature is available.  Must be non-zero. */
#undef BUILDOPT_OMP_FFTW
#define BUILDOPT_OMP_FFTW 1


/* Signalizes the feature is available.  Must be non-zero. */
#undef BUILDOPT_ISFINITE
#define BUILDOPT_ISFINITE 1


#undef BUILDOPT_PRECISION_SINGLE
#define BUILDOPT_PRECISION_SINGLE 1


#undef BUILDOPT_PRECISION_DOUBLE
#define BUILDOPT_PRECISION_DOUBLE 2


#undef BUILDOPT_PRECISION_QUAD
#define BUILDOPT_PRECISION_QUAD 3


#undef BUILDOPT_MPI
#define BUILDOPT_MPI 1


#undef BUILDOPT_SIMD_NONE
#define BUILDOPT_SIMD_NONE 0


#undef BUILDOPT_SIMD_AVX
#define BUILDOPT_SIMD_AVX 1


#undef BUILDOPT_SIMD_AVX2
#define BUILDOPT_SIMD_AVX2 2


#undef BUILDOPT_SIMD_AVX512
#define BUILDOPT_SIMD_AVX512 3


#undef BUILDOPT_SIMD_NEON
#define BUILDOPT_SIMD_NEON 4


#undef BUILDOPT_SIMD_SSE41
#define BUILDOPT_SIMD_SSE41 5


#undef LIB_NA_STR
#define LIB_NA_STR "n/a"


#undef LIB_NA_VAL
#define LIB_NA_VAL (-1)


#ifdef __cplusplus
}
#endif


#endif
