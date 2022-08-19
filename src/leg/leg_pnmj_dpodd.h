/* This header file is not a part of API. */


#ifndef __LEG_PNMJ_DPODD_H__
#define __LEG_PNMJ_DPODD_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(leg_pnmj_dpodd)(unsigned long, const REAL *, REAL *,
                                  REAL *, const int *, int *, int *,
                                  CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
