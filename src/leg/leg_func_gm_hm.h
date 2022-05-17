/* This header file is not a part of API. */


#ifndef __LEG_FUNC_GM_HM_H__
#define __LEG_FUNC_GM_HM_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** Computes the ``gm`` and ``hm`` coefficients that are useful to evaluate the
 * integrals of fully-normalized associated Legendre functions.  The
 * coefficients are defined as (see Eqs. A.2 and A.12 of Fukushima, 2014)
 *
 * \verbatim embed:rst:leading-asterisk
 *  .. math::
 *
 *      g_m &= \frac{1}{2 (m + 1)} \sqrt{\frac{m (2m + 1) (2m - 1)}{m - 1}}{,} 
 *      \\
 *      h_m &= \frac{m - 2}{m + 1}{,}
 *
 * \endverbatim
 *
 * for all ``m = 0``, ``1``, ``...``, ``nmax``.
 *
 *
 * \verbatim embed:rst:leading-asterisk
 *      .. note::
 *
 *          All coefficients ``gm`` for ``n < 3`` are set to zero, provided 
 *          that they do not exceed the maximum harmonic degree ``nmax``.
 *
 * \endverbatim
 *
 * **References**:
 *
 * * Fukushima T (2014) Numerical computation of spherical harmonics of 
 *   arbitrary degree and order by extending exponent of floating point 
 *   numbers.  Computers and Geosciences 63, 17-21, doi: 
 *   10.1016/j.cageo.2013.10.010
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
 * @param[out] gm Output coefficients, the coefficient \f$g_m\f$ is stored as 
 * ``gm[m]``.  For a given ``nmax``, there is ``nmax + 1`` elements of ``gm``.
 *
 * @param[out] hm Output coefficients, the coefficient \f$h_m\f$ is stored as 
 * ``hm[m]``.  For a given ``nmax``, there is ``nmax + 1`` elements of ``hm``.
 *
 * */
extern void CHARM(leg_func_gm_hm)(unsigned long nmax,
                                  const REAL *r,
                                  const REAL *ri,
                                  REAL *gm,
                                  REAL *hm);


#ifdef __cplusplus
}
#endif


#endif
