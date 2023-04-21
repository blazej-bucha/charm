/* This header file is not a part of API. */


#ifndef __CHECK_SHC_DDAV_H__
#define __CHECK_SHC_DDAV_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int check_shc_ddav(void (*)(const CHARM(shc) *,
                                        const CHARM(shc) *,
                                        unsigned long,
                                        REAL *,
                                        CHARM(err) *));


#ifdef __cplusplus
}
#endif


#endif
