/* This header file is not a part of API. */


#ifndef __SHS_GRD_FFT_LC_H__
#define __SHS_GRD_FFT_LC_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "shs_lc_struct.h"
#include "shs_max_npar.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_grd_fft_lc)(unsigned long,
                                  REAL,
                                  int,
                                  CHARM(lc) *,
                                  _Bool,
                                  REAL_SIMD *,
                                  int,
                                  size_t,
                                  REAL *,
                                  REAL *);


#ifdef __cplusplus
}
#endif


#endif
