/* This header file is not a part of API. */


#ifndef __MISC_STR2REAL_H__
#define __MISC_STR2REAL_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern REAL CHARM(misc_str2real)(const char *,
                                 const char *,
                                 CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
