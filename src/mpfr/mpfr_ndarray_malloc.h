/* This header file is not a part of API. */


#ifndef __MPFR_NDARRAY_MALLOC_H__
#define __MPFR_NDARRAY_MALLOC_H__


#include <config.h>
#include "../prec.h"
#include <mpfr.h>
#include "mpfr_ndarray.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern mpfr_ndarray *CHARM(mpfr_ndarray_malloc)(mpfr_prec_t,
                                                size_t,
                                                ...);


#ifdef __cplusplus
}
#endif


#endif
