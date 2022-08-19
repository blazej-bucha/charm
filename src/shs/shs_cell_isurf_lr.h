/* This header file is not a part of API. */


#ifndef __SHS_CELL_ISURF_LR_H__
#define __SHS_CELL_ISURF_LR_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_cell_isurf_lr)(REAL, REAL, size_t,
                                     REAL_SIMD, REAL_SIMD,
                                     REAL_SIMD, REAL_SIMD,
                                     unsigned long, unsigned long,
                                     REAL *);


#ifdef __cplusplus
}
#endif


#endif
