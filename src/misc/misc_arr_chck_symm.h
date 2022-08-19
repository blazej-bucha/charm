/* This header file is not a part of API. */


#ifndef __MISC_ARR_CHCK_SYMM_H__
#define __MISC_ARR_CHCK_SYMM_H__

#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/**
 * Checks whether a double ``arr`` array of ``size`` is symmetric with respect
 * to ``center``, that is, meets the following criteria:
 *
 * 1) ``arr[i] - center == center - arr[size - 1 - i]`` for ``i = 0, 1, ..., 
 * floor(size / 2)`` (up to some numerical threshold ``eps``),
 *
 * 2) if ``size`` is odd, the middle element of ``arr`` is equal to ``center`` 
 * (within some numerical threshold ``eps``).
 *
 * @return ``0`` if all conditions are met; otherwise, returned is the number 
 * of condition that is not satisfied.  In case of an error, returned is the 
 * ``-9999`` value.
 *
 * Error reported by the function (if any) is written to ``err``.
 *
 * */
extern int CHARM(misc_arr_chck_symm)(const REAL *arr, size_t size, REAL center,
                                     REAL eps, CHARM(err) *err);


#ifdef __cplusplus
}
#endif


#endif
