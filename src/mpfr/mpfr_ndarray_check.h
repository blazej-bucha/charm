/* This header file is not a part of API. */


#ifndef __MPFR_NDARRAY_CHECK_H__
#define __MPFR_NDARRAY_CHECK_H__


#include <config.h>
#include "../prec.h"
#include <mpfr.h>
#include "mpfr_ndarray.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern int CHARM(mpfr_ndarray_check)(const mpfr_ndarray *,
                                     size_t,
                                     ...);


#ifdef __cplusplus
}
#endif


#endif
