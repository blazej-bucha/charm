/* This header file is not a part of API. */


#ifndef __SHS_CELL_KERNEL_H__
#define __SHS_CELL_KERNEL_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


void CHARM(shs_cell_kernel)(unsigned long, unsigned long,
                            const CHARM(shc) *,
                            const REAL *, const REAL *,
                            REAL, REAL,
                            REAL, REAL,
                            REAL, REAL,
                            const REAL *, const REAL *,
                            const int *, const int *,
                            REAL *, REAL *, REAL *,
                            const REAL *, const REAL *, const REAL *,
                            const REAL *, const REAL *, const REAL *,
                            const REAL *,
                            _Bool,
                            REAL *, REAL *,
                            REAL *, REAL *);


#ifdef __cplusplus
}
#endif


#endif
