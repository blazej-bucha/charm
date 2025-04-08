/* This header file is not a part of API. */


#ifndef __MPFR_FACT_H__
#define __MPFR_FACT_H__


#include <config.h>
#include "../prec.h"
#include <mpfr.h>


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(mpfr_fact)(unsigned x,
                             mpfr_t out);


#ifdef __cplusplus
}
#endif


#endif
