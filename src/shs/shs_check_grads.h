/* This header file is not a part of API. */


#ifndef __SHS_CHECK_GRADS_H__
#define __SHS_CHECK_GRADS_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_check_grads)(int,
                                   int,
                                   int,
                                   CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
