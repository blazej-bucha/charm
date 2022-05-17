/* This header file is not a part of API. */


#ifndef __LEG_FUNC_ANM_BNM_H__
#define __LEG_FUNC_ANM_BNM_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/* For a given harmonic order ``m <= nmax`` and all degrees ``n = m + 2``, ``m
 * + 3``, ``...``, ``nmax``, computes the ``anm`` and ``bnm`` coefficients of
 * the three-term recursion formula for fully-normalized associated Legendre
 * functions (see Eqs. 5 - 7 of Fukushima, 2012).
 *
 * \verbatim embed:rst:leading-asterisk
 *      .. note::
 *
 *          All coefficients for ``n = m`` and ``n = m + 1`` are set to zero, 
 *          provided that they do not exceed the maximum harmonic degree 
 *          ``nmax``.
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
 * @param[out] anm Coefficients of the recursion formula, the coefficient
 * \f$a_{nm}\f$ is stored as ``anm[n]``.  For a given ``nmax``, there is
 * ``nmax + 1`` elements of ``anm``.
 *
 * @param[out] bnm Coefficients of the recursion formula, the coefficient
 * \f$b_{nm}\f$ is stored as ``bnm[n]``.  For a given ``nmax``, there is
 * ``nmax + 1`` elements of ``bnm``.
 *
 * */
extern void CHARM(leg_func_anm_bnm)(unsigned long nmax,
                                    unsigned long m,
                                    const REAL *r,
                                    const REAL *ri,
                                    REAL *anm,
                                    REAL *bnm);


#ifdef __cplusplus
}
#endif


#endif
