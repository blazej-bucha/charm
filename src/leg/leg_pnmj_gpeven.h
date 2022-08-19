/* This header file is not a part of API. */


#ifndef __LEG_PNMJ_GPEVEN_H__
#define __LEG_PNMJ_GPEVEN_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(leg_pnmj_gpeven)(unsigned long, unsigned long, unsigned long,
                                   const REAL *, const REAL *, REAL *,
                                   const int *, const int *, int *,
                                   CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
