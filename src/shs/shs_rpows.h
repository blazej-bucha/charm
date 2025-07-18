/* This header file is not a part of API. */


#ifndef __SHS_RPOWS_H__
#define __SHS_RPOWS_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_rpows)(REAL_SIMD,
                             REAL_SIMD,
                             size_t,
                             size_t,
                             REAL_SIMD *);


#ifdef __cplusplus
}
#endif


#endif
