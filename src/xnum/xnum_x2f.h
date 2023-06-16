/* This header file is not a part of API. */


#ifndef __XNUM_X2F_H__
#define __XNUM_X2F_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** @brief Translates an X-number into an F-number.
 *
 * @details The function is due to Fukushima (2012), Table 6.
 *
 * @param[in] x  Significand of the X-number
 * @param[in] ix Exponent of the X-number
 *
 * @returns The F-number representation of an input X-number.
 *
 * */
extern REAL CHARM(xnum_x2f)(REAL x,
                            int ix);


#ifdef __cplusplus
}
#endif


#endif
