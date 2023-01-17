/* This header file is not a part of API. */


#ifndef __MISC_POLAR_OPTIMIZATION_APPLY_H__
#define __MISC_POLAR_OPTIMIZATION_APPLY_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** @param[in] m Spherical harmonic order.
 *
 * @param[in] nmax Maximum spherical harmonic degree of the synthesis/analysis.
 *
 * @param[in] sinlat Sines of the latitudes.
 *
 * @param[in] threshold Polar optimization threshold.
 *
 * @return True if polar optimization can be applied, False otherwise.
 *
 * */
extern _Bool CHARM(misc_polar_optimization_apply)(unsigned long m,
                                                  unsigned long nmax,
                                                  REAL_SIMD sinlat,
                                                  REAL_SIMD threshold);


#ifdef __cplusplus
}
#endif


#endif
