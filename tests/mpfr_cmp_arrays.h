/* This header file is not a part of API. */


#ifndef __MPFR_CMP_ARRAYS_H__
#define __MPFR_CMP_ARRAYS_H__


#include <config.h>
#include <mpfr.h>


#ifdef __cplusplus
extern "C"
{
#endif


extern long int mpfr_cmp_arrays(mpfr_t *,
                                mpfr_t *,
                                size_t,
                                mpfr_t);


#ifdef __cplusplus
}
#endif


#endif
