/* This header file is not a part of API. */


#ifndef __LEG_FUNC_R_RI_H__
#define __LEG_FUNC_R_RI_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/** Computes auxiliary constants ``r`` and ``ri`` given as
 *
 * \verbatim embed:rst:leading-asterisk
 *  .. math::
 *
 *      r_m           &= \sqrt{m}{,}\\
 *      \mathrm{ri}_m &= \dfrac{1}{r_m}{,}
 *
 * \endverbatim
 *
 * for all ``m = 0``, ``1``, ``...``, ``2 * nmax + 3``.
 *
 * \verbatim embed:rst:leading-asterisk
 *      .. note::
 *
 *          Coefficients for ``m = 0`` are set to zero.
 *
 * \endverbatim
 *
 * \verbatim embed:rst:leading-asterisk
 *      .. warning::
 *
 *          For a given ``nmax``, the arrays ``r`` and ``ri`` must have access 
 *          to ``2 * nmax + 4`` elements.
 *
 * \endverbatim
 *
 * This function is useful to prepare the inputs to
 * ``CHARM(leg_func_anm_bnm)``, ``CHARM(leg_func_dm)`` and
 * ``CHARM(leg_func_gm_hm)``.
 *
 * @param[in] nmax Maximum harmonic degree.
 *
 * @param[out] r The coefficient \f$r_m\f$ can be accessed as ``r[m]``.  For a 
 * given ``nmax``, there is ``2 * nmax + 4`` elements of ``r``.
 *
 * @param[out] ri The coefficient \f$\mathrm{ri}_m\f$ can be accessed as 
 * ``ri[m]``.  For a given ``nmax``, there is ``2 * nmax + 4`` elements of 
 * ``ri``.
 *
 * */
extern void CHARM(leg_func_r_ri)(unsigned long nmax,
                                 REAL *r,
                                 REAL *ri);


#ifdef __cplusplus
}
#endif


#endif
