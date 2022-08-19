/* This header file is not a part of API. */


#ifndef __MISC_ARR_CHCK_LIN_INCR_H__
#define __MISC_ARR_CHCK_LIN_INCR_H__

#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/**
 * Checks whether elements of a double ``a`` array of ``size`` meet the
 * following conditions:
 *
 * 1) the ``first`` element of the array is smaller than or equal to the
 *    smallest element from the elements that follow it,
 *
 * 2) the difference ``a[i + every_nth] - a[i]`` for ``i = first``,
 *    ``first + every_nth``, ..., ``size - 1`` is positive and constant (within
 *    some numerical threshold ``eps``).
 *
 * For instance, if ``first = 0`` and ``every_nth = 1``, the function checks
 * whether all elements of ``a`` are linearly increasing with a constant step.
 * If ``every_nth = 2``, the function does the same, but for ``every_nth``
 * element starting from ``first``.
 *
 * @return ``0`` if all conditions are met; otherwise, returned is the number 
 * of the condition that is not satisfied.  In case of an error, returned is 
 * the ``-9999`` value.
 *
 * Error reported by the function (if any) is written to ``err``.
 *
 * */
extern int CHARM(misc_arr_chck_lin_incr)(const REAL *a, size_t na,
                                         size_t first, size_t every_nth,
                                         REAL eps, CHARM(err) *err);


#ifdef __cplusplus
}
#endif


#endif
