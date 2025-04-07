/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#if CHARM_QUAD
#   define MPFR_WANT_FLOAT128
#   define mpfr_float128 __float128
#endif
#include "../mpfr/mpfr_get_real.h"
#include "../mpfr/mpfr_set_real.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../mpfr/mpfr_check_bits.h"
/* ------------------------------------------------------------------------- */






/* This function is not a part of the API and is used with PyHarm only.
 *
 * See "gfm_cap_q_check_prec_pywrap.c" for details. */
CHARM_EXTERN void CHARM_CDECL CHARM(gfm_cap_q_pywrap)(REAL rref,
                                                      REAL r,
                                                      REAL psi,
                                                      unsigned long nmax,
                                                      unsigned pmax,
                                                      unsigned kmin,
                                                      unsigned kmax,
                                                      unsigned imax,
                                                      int zone,
                                                      unsigned type,
                                                      int NBITS,
                                                      REAL *qkpin,
                                                      CHARM(err) *err)
{
    /* Checks */
    /* --------------------------------------------------------------------- */
    mpfr_prec_t NBITS_mpfr = NBITS;
    CHARM(mpfr_check_bits)(NBITS_mpfr, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* --------------------------------------------------------------------- */


    /* Compute the truncation coefficients */
    /* --------------------------------------------------------------------- */
    mpfr_t rref_mpfr, r_mpfr, psi_mpfr;
    mpfr_inits2(NBITS_mpfr, rref_mpfr, r_mpfr, psi_mpfr, (mpfr_ptr)NULL);


    mpfr_set_REAL(rref_mpfr, rref, MPFR_RNDN);
    mpfr_set_REAL(r_mpfr, r, MPFR_RNDN);
    mpfr_set_REAL(psi_mpfr, psi, MPFR_RNDN);


    mpfr_t *qkpin_mpfr = NULL;


    size_t size = CHARM(gfm_cap_nq)(nmax, pmax, kmin, kmax, imax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }


    qkpin_mpfr = (mpfr_t *)malloc(size * sizeof(mpfr_t));
    if (qkpin_mpfr == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (size_t v = 0; v < size; v++)
        mpfr_init2(qkpin_mpfr[v], NBITS);


    CHARM(gfm_cap_q)(rref_mpfr, r_mpfr, psi_mpfr,
                     nmax, pmax, kmin, kmax, imax,
                     zone, type, NBITS_mpfr,
                     qkpin_mpfr, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */


    /* Type cast the truncation coefficients to "REAL" */
    /* --------------------------------------------------------------------- */
    for (size_t v = 0; v < size; v++)
        qkpin[v] = mpfr_get_REAL(qkpin_mpfr[v], MPFR_RNDN);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
EXIT:
    if (qkpin_mpfr != NULL)
        for (size_t v = 0; v < size; v++)
            mpfr_clear(qkpin_mpfr[v]);
    free(qkpin_mpfr);


    mpfr_clears(rref_mpfr, r_mpfr, psi_mpfr, (mpfr_ptr)NULL);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
    /* --------------------------------------------------------------------- */
}

