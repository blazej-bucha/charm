/* This header file is not a part of API. */


#ifndef __CHECK_CRD_POINT_QUAD_H__
#define __CHECK_CRD_POINT_QUAD_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int check_crd_point_quad(CHARM(point) *(*)(unsigned long,
                                                       REAL));


#ifdef __cplusplus
}
#endif


#endif
