/* This header file is not a part of API. */


#ifndef __SHS_GRD_LR_H__
#define __SHS_GRD_LR_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "shs_lc_struct.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_grd_lr)(unsigned long,
                              REAL,
                              REAL,
                              size_t,
                              int,
                              int,
                              size_t,
                              CHARM(lc) *,
                              _Bool,
                              REAL *,
                              REAL *);


#ifdef __cplusplus
}
#endif


#endif
