/* This header file is not a part of API. */


#ifndef __SHC_WRITE_MTDT_H__
#define __SHC_WRITE_MTDT_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


void CHARM(shc_write_mtdt)(unsigned long, REAL, REAL, const char *, FILE *,
                           CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
