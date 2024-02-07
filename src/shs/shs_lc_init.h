/* This header file is not a part of API. */


#ifndef __SHS_LC_INIT_H__
#define __SHS_LC_INIT_H__


#include "../prec.h"
#include "../simd/simd.h"


/* This routine *must* be called after the "CHARM(lc)" structure is declared
 * and before it is used for the first time. */
void CHARM(shs_lc_init)(CHARM(lc) *);


#endif
