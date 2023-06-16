/* This header file is not a part of API. */


#ifndef __LEG_PNMJ_LENGTH_H__
#define __LEG_PNMJ_LENGTH_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** Returns the total number of Fourier coefficients (to be) initialized by the
 * ``charm@P@_leg_pnmj_calloc()`` and ``charm@P@_leg_pnmj_malloc()`` functions
 * for a maximum harmonic degree ``nmax``.
 *
 * */
extern size_t CHARM(leg_pnmj_length)(unsigned long nmax);



#ifdef __cplusplus
}
#endif


#endif
