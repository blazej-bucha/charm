/* This header file is not a part of API. */


#ifndef __SHS_POINT_KERNEL_H__
#define __SHS_POINT_KERNEL_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


void CHARM(shs_point_kernel)(unsigned long, unsigned long, const CHARM(shc) *,
                             const REAL *, const REAL *, REAL,
                             const REAL *, const int *, const REAL *,
                             const REAL *, _Bool,
                             REAL *, REAL *, REAL *, REAL *);


#ifdef __cplusplus
}
#endif


#endif
