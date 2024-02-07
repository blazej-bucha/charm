/* This header file is not a part of API. */


#ifndef __LEG_FUNC_ENM_H__
#define __LEG_FUNC_ENM_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/* For a given harmonic order ``m <= nmax`` and all degrees ``n = m``, ``m
 * + 1``, ``...``, ``nmax``, computes the ``enm`` coefficients of the singular
 * three-term fix-order recursion formula for the first-order derivatives of
 * the fully-normalized associated Legendre functions (see Eqs. 19 of
 * Fukushima, 2012).
 *
 * \verbatim embed:rst:leading-asterisk
 *      .. note::
 *
 *          All coefficients for ``n = m`` are set to zero, provided that they
 *          do not exceed the maximum harmonic degree ``nmax``.
 *
 * \endverbatim
 *
 * **References**:
 *
 * * Fukushima T (2012) Numerical computation of spherical harmonics of 
 *   arbitrary degree and order by extending exponent of floating point
 *   numbers: II first-, second-, and third-order derivatives.  Doi:
 *   10.1007/s00190-012-0561-8
 *
 * @param[in] nmax Maximum harmonic degree of the coefficients.
 *
 * @param[in] m Harmonic order.
 *
 * @param[in] r Some auxiliary coefficients pre-computed by
 * ``CHARM(leg_func_r_ri)``.  Such pre-computation is useful when calling the
 * function multiple times for various ``m``.
 *
 * @param[in] ri Some auxiliary coefficients pre-computed by
 * ``CHARM(leg_func_r_ri)``.  Such pre-computation is useful when calling the
 * function multiple times for various ``m``.
 *
 * @param[out] enm Coefficients of the recursion formula, the coefficient
 * \f$e_{nm}\f$ is stored as ``enm[n]``.  For a given ``nmax``, there is
 * ``nmax + 1`` elements of ``enm``.
 *
 * */
extern void CHARM(leg_func_enm)(unsigned long nmax,
                                unsigned long m,
                                const REAL *r,
                                const REAL *ri,
                                REAL *enm);


#ifdef __cplusplus
}
#endif


#endif
