/* This header file is not a part of API. */


#ifndef __SHC_ARITHMETICS_CHECKS_H__
#define __SHC_ARITHMETICS_CHECKS_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shc_arithmetics_checks)(CHARM(shc) *,
                                          CHARM(shc) *,
                                          CHARM(shc) *,
                                          unsigned long,
                                          unsigned long,
                                          CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
