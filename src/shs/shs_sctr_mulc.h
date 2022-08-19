/* This header file is not a part of API. */


#ifndef __SHS_SCTR_MULC_H__
#define __SHS_SCTR_MULC_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_sctr_mulc)(size_t, size_t, REAL_SIMD, REAL_SIMD,
                                 REAL *, REAL_SIMD, REAL *);


#ifdef __cplusplus
}
#endif


#endif

