/* This header file is not a part of API. */


#ifndef __XNUM_XLSUM2_H__
#define __XNUM_XLSUM2_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** @brief Computes a two-term linear sum of X-numbers with F-number
 * coefficients.
 *
 * @details The function is due to Fukushima (2012), Table 8.
 *
 *
 * @param[in] f   F-number coefficient of the first linear term
 *
 * @param[in] x   Significand of the first X-number term
 *
 * @param[in] g   F-number coefficient of the second linear term
 *
 * @param[in] y   Significand of the second X-number term
 *
 * @param[in] ix  Exponent of the first X-number term
 *
 * @param[in] iy  Exponent of the second X-number term
 *
 * @param[out] z  Pointer to the significand of the output X-number
 *
 * @param[out] iz Pointer to the exponent of the output X-number
 *
 * */
extern void CHARM(xnum_xlsum2)(REAL f,
                               REAL x,
                               REAL g,
                               REAL y,
                               REAL *z,
                               int ix,
                               int iy,
                               int *iz);


#ifdef __cplusplus
}
#endif


#endif
