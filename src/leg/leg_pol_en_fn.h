/* This header file is not a part of API. */


#ifndef __LEG_POL_EN_FN_H__
#define __LEG_POL_EN_FN_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** Computes the ``en`` and ``fn`` coefficients of the Bonnet's recursion
 * formula for un-normalized Legendre polynomials up to degree ``nmax``:
 *
 * \verbatim embed:rst:leading-asterisk
 *  .. math::
 *
 *      P_n(t) = e_n \, t \, P_{n - 1}(t) - f_n \, P_{n - 2}(t), \quad t \in 
 *      [-1, 1] {.}
 *
 * \endverbatim
 *
 * \verbatim embed:rst:leading-asterisk
 *      .. note::
 *
 *          All coefficients for ``n <= 1`` are set to zero, provided that they 
 *          do not exceed the maximum harmonic degree ``nmax``.
 *
 * \endverbatim
 *
 * @param[in] nmax Maximum harmonic degree of the coefficients.
 *
 * @param[out] en Coefficients of the recursion formula related to the \f$P_{n
 * - 1}\f$ term, the coefficient \f$e_n\f$ is stored as ``e[n]``. For a given
 *   ``nmax``, there is ``nmax + 1`` elements of ``e``.
 *
 * @param[out] fn Coefficients of the recursion formula related to the \f$P_{n 
 * - 2}\f$ term, the coefficient \f$f_n\f$ is stored as ``f[n]``.  For a given 
 *   ``nmax``, there is ``nmax + 1`` elements of ``f``.
 *
 * */
extern void CHARM(leg_pol_en_fn)(unsigned long nmax,
                                 REAL *en,
                                 REAL *fn);


#ifdef __cplusplus
}
#endif


#endif
