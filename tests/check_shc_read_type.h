/* This header file is not a part of API. */


#ifndef __CHECK_SHC_READ_TYPE_H__
#define __CHECK_SHC_READ_TYPE_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int check_shc_read_type(unsigned long (*)(const char *,
                                                      unsigned long,
                                                      CHARM(shc) *,
                                                      CHARM(err) *));


#ifdef __cplusplus
}
#endif


#endif
