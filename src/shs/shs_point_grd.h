/* This header file is not a part of API. */


#ifndef __SHS_POINT_GRD_H__
#define __SHS_POINT_GRD_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


void CHARM(shs_point_grd)(const CHARM(crd) *, const CHARM(shc) *,
                          unsigned long, REAL *, CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
