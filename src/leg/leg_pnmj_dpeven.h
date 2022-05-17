/* This header file is not a part of API. */


#ifndef __LEG_PNMJ_DPEVEN_H__
#define __LEG_PNMJ_DPEVEN_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


void CHARM(leg_pnmj_dpeven)(unsigned long, const REAL *, REAL *,
                            REAL *, const int *, int *, int *,
                            CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
