/* This header file is not a part of API. */


#ifndef __LEG_FUNC_DM_H__
#define __LEG_FUNC_DM_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/* Computes the ``dm`` coefficients of sectorial fully-normalized associated 
 * Legendre functions (see Eq. 9 of Fukushima, 2012) up to degree ``nmax``.
 *
 * \verbatim embed:rst:leading-asterisk
 *      .. note::
 *
 *          The coefficient for ``n = 0`` is set to zero.
 *
 * \endverbatim
 *
 * **References**:
 *
 * * Fukushima T (2012) Numerical computation of spherical harmonics of 
 *   arbitrary degree and order by extending exponent of floating point 
 *   numbers. Journal of Geodesy 86:271–285, doi: 10.1007/s00190-011-0519-2
 *
 * @param[in] nmax Maximum harmonic degree of the coefficients.
 *
 * @param[in] r Some auxiliary coefficients pre-computed by
 * ``CHARM(leg_func_r_ri)``.  Such pre-computation is useful when calling the
 * function multiple times for various ``m``.
 *
 * @param[in] ri Some auxiliary coefficients pre-computed by
 * ``CHARM(leg_func_r_ri)``.  Such pre-computation is useful when calling the
 * function multiple times for various ``m``.
 *
 * @param[out] dm Output coefficients, the coefficient \f$d_{m}\f$ is stored as 
 * ``dm[m]``.  For a given ``nmax``, there is ``nmax + 1`` elements of ``dm``.
 *
 * */
extern void CHARM(leg_func_dm)(unsigned long nmax,
                               const REAL *r,
                               const REAL *ri,
                               REAL *dm);


#ifdef __cplusplus
}
#endif


#endif
