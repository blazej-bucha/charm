/* This header file is not a part of API. */


#ifndef __SHS_CHECK_SINGLE_DERIVATIVE_H__
#define __SHS_CHECK_SINGLE_DERIVATIVE_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


#undef SHS_MAX_DERIVATIVE
#define SHS_MAX_DERIVATIVE 2


extern void CHARM(shs_check_single_derivative)(int,
                                               int,
                                               int,
                                               CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
