/* This header file is not a part of API. */


#ifndef __MPFR_BINOMIAL_H__
#define __MPFR_BINOMIAL_H__


#include <config.h>
#include "../prec.h"
#include <mpfr.h>
#include "mpfr_ndarray.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(mpfr_binomial)(mpfr_ndarray *,
                                 unsigned,
                                 mpfr_prec_t,
                                 CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
