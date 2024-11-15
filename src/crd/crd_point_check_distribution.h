/* This header file is not a part of API. */


#ifndef __CRD_POINT_CHECK_DISTRIBUTION_H__
#define __CRD_POINT_CHECK_DISTRIBUTION_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern _Bool CHARM(crd_point_check_distribution)(const CHARM(point) *,
                                                 CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
