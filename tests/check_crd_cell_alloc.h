/* This header file is not a part of API. */


#ifndef __CHECK_CRD_CELL_ALLOC_H__
#define __CHECK_CRD_CELL_ALLOC_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int check_crd_cell_alloc(CHARM(cell) *(*)(int,
                                                      size_t,
                                                      size_t));


#ifdef __cplusplus
}
#endif


#endif
