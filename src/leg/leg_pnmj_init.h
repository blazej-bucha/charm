/* This header file is not a part of API. */


#ifndef __LEG_PNMJ_INIT_H__
#define __LEG_PNMJ_INIT_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** Creates a "charm_pnmj" structure up to degree "nmax" using "ordering" and
 *  puts the "pnmj_coeffs" array to "charm_pnmj->pnmj".  */
extern CHARM(pnmj) *CHARM(leg_pnmj_init)(unsigned long nmax, 
                                         int ordering,
                                         REAL *pnmj_coeffs);


#ifdef __cplusplus
}
#endif


#endif
