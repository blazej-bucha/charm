/* This header file is not a part of API. */


#ifndef __MPFR_LEGENDRE_H__
#define __MPFR_LEGENDRE_H__


#include <config.h>
#include "../prec.h"
#include <mpfr.h>
#include "mpfr_ndarray.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(mpfr_legendre)(mpfr_ndarray *,
                                 unsigned long,
                                 unsigned long,
                                 const mpfr_t,
                                 mpfr_prec_t,
                                 CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
