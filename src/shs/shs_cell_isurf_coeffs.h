/* This header file is not a part of API. */


#ifndef __SHS_CELL_ISURF_COEFFS_H__
#define __SHS_CELL_ISURF_COEFFS_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


void CHARM(shs_cell_isurf_coeffs)(const CHARM(shc) *, unsigned long,
                                  const CHARM(shc) *, unsigned long,
                                  unsigned long, unsigned long,
                                  REAL *, REAL *, REAL *,
                                  REAL *, CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
