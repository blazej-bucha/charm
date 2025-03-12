/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if CHARM_QUAD
#   define MPFR_WANT_FLOAT128
#   define _Float128 __float128
#endif
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_propagate.h"
#include "../mpfr/mpfr_check_bits.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
/* ------------------------------------------------------------------------- */






/* This function is not a part of the API and is used with PyHarm only.
 *
 * The main rationale is to provide a Python wrapper for
 * "gfm_cap_q_check_prec".  "rref", "r" and "psi" are represented by "REAL" and
 * then are converted to "mpfr_t" (the conversion is likely not exact,
 * depending on the values of these variables).  To represent "mpfr_prec_t",
 * "NBITS" and "NBITSREF" are assumed to be "int".  Note that "mpfr_prec_t" can
 * also be "long" on some (perhaps most) architectures.  The main problem here
 * is that we do not know which representation of "mpfr_prec_t" is relevant for
 * the host system.  To be on the safe side, "int" is used.  As a result, this
 * wrapper doesn't allow to use rather extremely large numbers of bits.  If
 * this is what you need, use "gfm_cap_q_check_prec" from CHarm without any
 * limitations.  In any case, "int" should allow you to use sufficiently large
 * integers on most modern systems, so that most practical applications should
 * not be limited by this.  It seems MPFR never defines "mpfr_prec_t" as
 * "short" so, "int" should be OK.
 *
 * Using "int" to represent "NBITS" and "NBITSREF" is the main downside of this
 * wrapper (but hopefully is better than nothing for Python users).  Other
 * approaches are possible with external Python libraries, but we do not want
 * to add further dependencies.
 *
 * Since PyHarm does not support quadruple precision, the installation scripts
 * do not compile this function when quadruple precision is used.  However, the
 * function itself is written such that it could be compiled even in quadruple
 * precision.  This, however, requires to have "MPFR" installed with the
 * "--enable-float128" installation flag which is why this function is not
 * compiled if using quadruple precision. */
CHARM_EXTERN long CHARM_CDECL CHARM(gfm_cap_q_check_prec_pywrap)(
                                        REAL rref,
                                        REAL r,
                                        REAL psi,
                                        unsigned long nmax,
                                        unsigned pmax,
                                        unsigned kmin,
                                        unsigned kmax,
                                        unsigned imax,
                                        unsigned type,
                                        int NBITS,
                                        int NBITSREF,
                                        CHARM(err) *err)
{
    /* Checks */
    /* --------------------------------------------------------------------- */
    CHARM(mpfr_check_bits)(NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return -1;
    }


    CHARM(mpfr_check_bits)(NBITSREF, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return -2;
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    mpfr_prec_t NBITS_mpfr    = NBITS;
    mpfr_prec_t NBITSREF_mpfr = NBITSREF;


    mpfr_t rref_mpfr, r_mpfr, psi_mpfr;
    mpfr_inits2(NBITS_mpfr, rref_mpfr, r_mpfr, psi_mpfr, (mpfr_ptr)NULL);


#if CHARM_FLOAT
    mpfr_set_flt(rref_mpfr, rref, MPFR_RNDN);
    mpfr_set_flt(r_mpfr, r, MPFR_RNDN);
    mpfr_set_flt(psi_mpfr, psi, MPFR_RNDN);
#elif CHARM_QUAD
    mpfr_set_float128(rref_mpfr, rref, MPFR_RNDN);
    mpfr_set_float128(r_mpfr, r, MPFR_RNDN);
    mpfr_set_float128(psi_mpfr, psi, MPFR_RNDN);
#else
    mpfr_set_d(rref_mpfr, rref, MPFR_RNDN);
    mpfr_set_d(r_mpfr, r, MPFR_RNDN);
    mpfr_set_d(psi_mpfr, psi, MPFR_RNDN);
#endif


    long ret = -3;
    ret = CHARM(gfm_cap_q_check_prec)(rref_mpfr,
                                      r_mpfr,
                                      psi_mpfr,
                                      nmax,
                                      pmax,
                                      kmin,
                                      kmax,
                                      imax,
                                      type,
                                      NBITS_mpfr,
                                      NBITSREF_mpfr,
                                      err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    mpfr_clears(rref_mpfr, r_mpfr, psi_mpfr, (mpfr_ptr)NULL);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return ret;
    /* --------------------------------------------------------------------- */
}

