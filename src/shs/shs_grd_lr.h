/* This header file is not a part of API. */


#ifndef __SHS_GRD_LR_H__
#define __SHS_GRD_LR_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_grd_lr)(unsigned long, REAL, REAL, size_t, int,
                              REAL_SIMD, REAL_SIMD, REAL_SIMD, REAL_SIMD,
                              _Bool, REAL *, REAL *);


#ifdef __cplusplus
}
#endif


#endif
