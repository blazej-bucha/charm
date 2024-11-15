/* This header file is not a part of API. */


#ifndef __SHS_RATIOS_H__
#define __SHS_RATIOS_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_ratios)(const REAL_SIMD *,
                              const REAL_SIMD *,
                              _Bool,
                              unsigned long,
                              REAL_SIMD *,
                              REAL_SIMD *);


#ifdef __cplusplus
}
#endif


#endif
