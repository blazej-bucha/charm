/* This header file is not a part of API. */


#ifndef __SHS_POINT_KERNEL_H__
#define __SHS_POINT_KERNEL_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_point_kernel)(unsigned long, unsigned long,
                                    const CHARM(shc) *, _Bool,
                                    const REAL *, const REAL *,
                                    REAL_SIMD *, const REAL *,
                                    const int *,
                                    REAL_SIMD *, REAL_SIMD *,
                                    REAL_SIMD *, REAL_SIMD *,
                                    REAL_SIMD *, REAL_SIMD *,
                                    REAL_SIMD *, REAL_SIMD *, REAL_SIMD *);


#ifdef __cplusplus
}
#endif


#endif
