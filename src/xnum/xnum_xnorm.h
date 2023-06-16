/* This header file is not a part of API. */


#ifndef __XNUM_XNORM_H__
#define __XNUM_XNORM_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** @brief Weakly normalizes an X-number.
 *
 * @details The function is due to Fukushima (2012), Table 7.
 *
 * @param[in]   x Pointer to the significand of the input X-number
 *
 * @param[in]  ix Pointer to the exponent of the input X-number
 *
 * @param[out]  x Pointer to the significand of the output X-number
 *
 * @param[out] ix Pointer to the exponent of the output X-number
 *
 * */
extern void CHARM(xnum_xnorm)(REAL *x,
                              int *ix);


#ifdef __cplusplus
}
#endif


#endif
