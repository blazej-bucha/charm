/* This header file is not a part of API. */


#ifndef __SHS_GET_MUR_DORDER_H__
#define __SHS_GET_MUR_DORDER_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_get_mur_dorder_npar)(const CHARM(shc) *shcs,
                                           int dr,
                                           int dlat,
                                           int dlon,
                                           REAL *mur,
                                           unsigned *dorder,
                                           size_t *npar,
                                           CHARM(err) *err);


#ifdef __cplusplus
}
#endif


#endif
