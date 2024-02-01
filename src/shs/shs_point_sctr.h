/* This header file is not a part of API. */


#ifndef __SHS_POINT_SCTR_H__
#define __SHS_POINT_SCTR_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_point_sctr)(const CHARM(point) *,
                                  const CHARM(shc) *,
                                  unsigned long,
                                  int,
                                  int,
                                  int,
                                  REAL **,
                                  CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
