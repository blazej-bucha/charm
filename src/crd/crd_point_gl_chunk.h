/* This header file is not a part of API. */


#ifndef __CRD_POINT_GL_CHUNK_H__
#define __CRD_POINT_GL_CHUNK_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern CHARM(point) *CHARM(crd_point_gl_chunk)(unsigned long,
                                               REAL,
                                               size_t,
                                               size_t,
                                               CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif

