/* This header file is not a part of API. */


#ifndef __MISC_POLAR_OPTIMIZATION_THRESHOLD_H__
#define __MISC_POLAR_OPTIMIZATION_THRESHOLD_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** @param[in] nmax Maximum degree of the synthesis/analysis.
 *
 * @return Threshold value to decide whether the polar optimization is to be
 * applied.
 *
 * */
extern REAL_SIMD CHARM(misc_polar_optimization_threshold)(unsigned long nmax);


#ifdef __cplusplus
}
#endif


#endif
