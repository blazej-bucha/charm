/* This header file is not a part of API. */


#ifndef __SHS_CELL_KERNEL_H__
#define __SHS_CELL_KERNEL_H__


#include <config.h>
#include <stdint.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_cell_kernel)(unsigned long, unsigned long,
                                   const CHARM(shc) *,
                                   const REAL *,
                                   const REAL *,
                                   REAL_SIMD,
                                   REAL_SIMD,
                                   REAL_SIMD,
                                   REAL_SIMD,
                                   REAL_SIMD,
                                   REAL_SIMD,
                                   const REAL *,
                                   const REAL *,
                                   const int64_t *,
                                   const int64_t *,
                                   REAL_SIMD *,
                                   REAL_SIMD *,
                                   REAL_SIMD *,
                                   const REAL *,
                                   const REAL *,
                                   const REAL *,
                                   const REAL *,
                                   const REAL *,
                                   REAL_SIMD,
                                   REAL_SIMD,
                                   REAL_SIMD,
                                   REAL_SIMD,
                                   REAL_SIMD,
                                   REAL_SIMD *,
                                   REAL_SIMD *,
                                   REAL_SIMD *,
                                   REAL_SIMD *);


#ifdef __cplusplus
}
#endif


#endif
