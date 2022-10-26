/* This header file is not a part of API. */


#ifndef __LEG_PNMJ_ALLOC_H__
#define __LEG_PNMJ_ALLOC_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern CHARM(pnmj) *CHARM(leg_pnmj_alloc)(unsigned long, int,
                                          void *(*)(size_t));


#ifdef __cplusplus
}
#endif


#endif
