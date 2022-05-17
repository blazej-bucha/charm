/* This header file is not a part of API. */


#ifndef __INTEG_SC_H__
#define __INTEG_SC_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/**
 * Analytically computes integrals:
 *
 * \verbatim embed:rst:leading-asterisk
 *      .. math::
 *
 *          s_j = \int_{u_0 + j \, \Delta u}^{u_0 + (j + 1) \, \Delta
 *          u} \sin(a_1 \, u) \, \cos(a_2 \, u) \, \mathrm{d} u
 *
 * \endverbatim
 *
 * for all integers ``0 <= j < nu`` by applying recurrences.
 *
 *
 * @param[in] u0 Lower integration limit for ``j = 0``.
 * @param[in] du Integration width.
 * @param[in] nu Number of integrals to be computed recursively, with 
 *   integration domain each time shifted by the integration width ``du``, ``j 
 *   = 0``, ``1``, ..., ``nu - 1``.
 * @param[in] a1 Real number.
 * @param[in] a2 Real number.
 * 
 * @param[out] s Pointer to an array of size ``nu`` with integrals ``s[j] = 
 * s_j``.
 *
 * */
extern void CHARM(integ_sc)(REAL u0,
                            REAL du,
                            size_t nu,
                            REAL a1,
                            REAL a2,
                            REAL *s);


#ifdef __cplusplus
}
#endif


#endif
