/* This header file is not a part of API. */


#ifndef __SHC_READ_MTDT_H__
#define __SHC_READ_MTDT_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shc_read_mtdt)(FILE *, unsigned long *, REAL *, REAL *,
                                 CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
