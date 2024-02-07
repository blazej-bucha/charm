/* This header file is not a part of API. */


#ifndef __SHS_CELL_CHECK_GRD_LONS_H__
#define __SHS_CELL_CHECK_GRD_LONS_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shs_cell_check_grd_lons)(const CHARM(cell) *,
                                           REAL *deltalon,
                                           CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
