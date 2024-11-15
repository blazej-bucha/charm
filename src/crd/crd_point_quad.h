/* This header file is not a part of API. */


#ifndef __CRD_POINT_QUAD_H__
#define __CRD_POINT_QUAD_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern CHARM(point) *CHARM(crd_point_quad)(unsigned long,
                                           REAL,
                                           void (*)(unsigned long,
                                                    size_t *,
                                                    size_t *),
                                           CHARM(point) *(*)(unsigned long,
                                                             REAL,
                                                             size_t,
                                                             size_t,
                                                             CHARM(err) *));


#ifdef __cplusplus
}
#endif


#endif

