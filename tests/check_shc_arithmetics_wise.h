/* This header file is not a part of API. */


#ifndef __CHECK_SHC_ARITHMETICS_WISE_H__
#define __CHECK_SHC_ARITHMETICS_WISE_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int check_shc_arithmetics_wise(void (*)(CHARM(shc) *,
                                                    const REAL *,
                                                    unsigned long,
                                                    unsigned long,
                                                    CHARM(err) *));


#ifdef __cplusplus
}
#endif


#endif
