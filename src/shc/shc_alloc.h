/* This header file is not a part of API. */


#ifndef __SHC_ALLOC_H__
#define __SHC_ALLOC_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern CHARM(shc) *CHARM(shc_alloc)(unsigned long, REAL, REAL,
                                    void *(*)(size_t));


#ifdef __cplusplus
}
#endif


#endif
