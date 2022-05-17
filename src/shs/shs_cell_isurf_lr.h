/* This header file is not a part of API. */


#ifndef __SHS_CELL_ISURF_LR_H__
#define __SHS_CELL_ISURF_LR_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


void CHARM(shs_cell_isurf_lr)(REAL, REAL, size_t,
                              REAL, REAL, REAL, REAL,
                              unsigned long, unsigned long,
                              REAL *);


#ifdef __cplusplus
}
#endif


#endif
