/* This header file is not a part of API. */


#ifndef __CRD_CELL_ALLOC_H__
#define __CRD_CELL_ALLOC_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern CHARM(cell) *CHARM(crd_cell_alloc)(int, size_t, size_t,
                                          void *(*)(size_t));


#ifdef __cplusplus
}
#endif


#endif

