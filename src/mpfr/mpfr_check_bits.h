/* This header file is not a part of API. */


#ifndef __MPFR_CHECK_BITS_H__
#define __MPFR_CHECK_BITS_H__


#include <config.h>
#include <mpfr.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(mpfr_check_bits)(mpfr_prec_t,
                                   CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
