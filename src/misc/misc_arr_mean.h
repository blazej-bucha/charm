/* This header file is not a part of API. */


#ifndef __MISC_ARR_MEAN_H__
#define __MISC_ARR_MEAN_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/**
 * @param[in] a A ``double`` array.
 *
 * @param[in] na Length of ``a``.
 *
 * @param[out] err Error reported by the function (if any).
 *
 * @return The mean value of ``a``.  In case of an error, returned is the 
 * ``nan`` value and the error is written to ``err``.
 *
 * */
extern REAL CHARM(misc_arr_mean)(const REAL *a, size_t na, CHARM(err) *err);


#ifdef __cplusplus
}
#endif


#endif
