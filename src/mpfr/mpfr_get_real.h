/* This header file is not a part of API. */


#ifndef __MPFR_GET_REAL_H__
#define __MPFR_GET_REAL_H__


#include <config.h>
#if CHARM_QUAD
#   define MPFR_WANT_FLOAT128
#endif
#include <mpfr.h>


/* Defines a macro to set a "mpfr_t" variable from a given "REAL" variable
 * based on the subroutines provided by MPFR. */
#if CHARM_FLOAT
#   define mpfr_get_REAL mpfr_get_flt
#elif CHARM_QUAD
#   define _Float128 __float128
#   define mpfr_get_REAL mpfr_get_float128
#else
#   define mpfr_get_REAL mpfr_get_d
#endif


#endif
