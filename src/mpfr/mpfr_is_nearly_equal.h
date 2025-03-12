/* This header file is not a part of API. */


#ifndef __MPFR_IS_NEARLY_EQUAL_H__
#define __MPFR_IS_NEARLY_EQUAL_H__


#include <config.h>
#include <mpfr.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern _Bool CHARM(mpfr_is_nearly_equal)(mpfr_t,
                                         mpfr_t,
                                         mpfr_t);


#ifdef __cplusplus
}
#endif


#endif
