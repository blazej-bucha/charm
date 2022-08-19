/* This header file is not a part of API. */


#ifndef __MISC_IS_NEARLY_EQUAL_H__
#define __MISC_IS_NEARLY_EQUAL_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/**
 * Judges whether two double numbers, ``a`` and ``b``, are equal up to some 
 * threshold ``eps``.
 *
 * **Example use**:
 *
 * \verbatim embed:rst:leading-asterisk
 * .. code-block:: c
 *
 *      if (charm@P@_misc_is_nearly_equal(a, b, eps))
 *      {
 *          // The two numbers are equal up to "eps"
 *      }
 *      else
 *      {
 *          // The two numbers are not equal up to "eps"
 *      }
 *
 * \endverbatim
 *
 * @return Boolean ``1`` if the two double numbers are equal up to the 
 *         threshold ``eps`` or ``0`` otherwise.
 *
 * */
extern _Bool CHARM(misc_is_nearly_equal)(REAL a, REAL b, REAL eps);


#ifdef __cplusplus
}
#endif


#endif
