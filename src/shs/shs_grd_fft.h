/* This header file is not a part of API. */


#ifndef __SHS_GRD_FFT_H__
#define __SHS_GRD_FFT_H__


#include <config.h>
#include <fftw3.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_grd_fft)(size_t, int, size_t, size_t,
                               REAL *, REAL *, REAL *, REAL,
                               FFTW(complex) *, FFTW(complex) *, size_t,
                               const REAL *, const REAL *,
                               REAL, const FFTW(plan), const REAL *,
                               REAL *, REAL *, REAL *);


#ifdef __cplusplus
}
#endif


#endif
