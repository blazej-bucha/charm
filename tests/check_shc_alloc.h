/* This header file is not a part of API. */


#ifndef __CHECK_SHC_ALLOC_H__
#define __CHECK_SHC_ALLOC_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int check_shc_alloc(CHARM(shc) *(*)(unsigned long,
                                                REAL,
                                                REAL));


#ifdef __cplusplus
}
#endif


#endif
