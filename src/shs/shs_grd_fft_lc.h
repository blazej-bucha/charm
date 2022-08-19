/* This header file is not a part of API. */


#ifndef __SHS_GRD_FFT_LC_H__
#define __SHS_GRD_FFT_LC_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_grd_fft_lc)(unsigned long, REAL, REAL_SIMD, REAL_SIMD,
                                  REAL_SIMD, REAL_SIMD, _Bool, REAL_SIMD, int,
                                  REAL *, REAL *);


#ifdef __cplusplus
}
#endif


#endif
