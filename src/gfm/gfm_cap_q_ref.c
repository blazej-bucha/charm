/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../mpfr/mpfr_ndarray.h"
#include "../mpfr/mpfr_ndarray_malloc.h"
#include "../mpfr/mpfr_ndarray_free.h"
#include "../mpfr/mpfr_fact.h"
#include "../mpfr/mpfr_check_bits.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "../misc/misc_idx_4d.h"
#include "gfm_check_p.h"
#include "gfm_check_kminkmax.h"
#include "gfm_cap_q_check_radius.h"
#include "gfm_cap_q_check_type.h"
#include "gfm_cap_q_ref.h"
/* ------------------------------------------------------------------------- */






/* Computes the independent reference values for the sum of near- and far-zone
 * truncation coefficients */
void CHARM(gfm_cap_q_ref)(const mpfr_t rref,
                          const mpfr_t r,
                          unsigned long nmax,
                          unsigned pmax,
                          unsigned kmin,
                          unsigned kmax,
                          unsigned imax,
                          unsigned type,
                          mpfr_prec_t NBITS,
                          mpfr_t *skpn,
                          CHARM(err) *err)
{
    /* Checks */
    /* --------------------------------------------------------------------- */
    CHARM(mpfr_check_bits)(NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_cap_q_check_radius)(rref, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_cap_q_check_radius)(r, NBITS, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_check_p)(pmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_check_kminkmax)(kmin, kmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(gfm_cap_q_check_type)(type, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    mpfr_ndarray *fact = CHARM(mpfr_ndarray_malloc)(NBITS, 1, pmax + 1);
    if (fact == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }
    for (unsigned p = 0; p <= pmax; p++)
        CHARM(mpfr_fact)(p, fact->data[p]);


    unsigned u = 0;
    unsigned v = 0;
    if (type == CHARM_GFM_Q00)
    {
        u = 0;
        v = 0;
    }
    else if (type == CHARM_GFM_Q10)
    {
        u = 1;
        v = 0;
    }
    else if (type == CHARM_GFM_Q11)
    {
        u = 1;
        v = 1;
    }
    else if (type == CHARM_GFM_Q20)
    {
        u = 2;
        v = 0;
    }
    else if (type == CHARM_GFM_Q21)
    {
        u = 2;
        v = 1;
    }
    else if (type == CHARM_GFM_Q22)
    {
        u = 2;
        v = 2;
    }


    size_t idx;
    mpfr_t r1k, rrefpow, t;
    mpfr_inits2(NBITS, r1k, rrefpow, t, (mpfr_ptr)NULL);


    mpfr_div(t, rref, r, MPFR_RNDN);
    for (unsigned k = kmin; k <= kmax; k++)
    {
        mpfr_pow_ui(rrefpow, rref, k + u, MPFR_RNDN);
        mpfr_si_div(r1k, ((k + u) % 2) ? -1 : 1, rrefpow, MPFR_RNDN);


        for (unsigned p = 1; p <= pmax; p++)
        {
            for (unsigned i = 0; i <= imax; i++)
            {
                idx = CHARM(misc_idx_4d)(k - kmin, p - 1, i, 0, pmax, imax + 1,
                                         nmax + 1);


#if CHARM_OPENMP
#pragma omp parallel default(shared)
#endif
                {
                    mpfr_t af, klpi, tpow, tmp1, tmp2;
                    mpfr_inits2(NBITS, af, klpi, tpow, tmp1, tmp2,
                                (mpfr_ptr)NULL);


                    {
                    unsigned long n;
#if CHARM_OPENMP
#pragma omp for
#endif
                    for (n = v; n <= nmax; n++)
                    {
                        mpfr_set_ui(af, 1, MPFR_RNDN);
                        for (unsigned e = v + 1; e <= k + u; e++)
                            mpfr_mul_ui(af, af, n + e, MPFR_RNDN);


                        mpfr_set_ui(klpi, 1, MPFR_RNDN);
                        mpfr_set_ui(tmp1, n + i + 4, MPFR_RNDN);
                        for (unsigned s = 2; s <= p; s++)
                        {
                            mpfr_sub_ui(tmp2, tmp1, s, MPFR_RNDN);
                            mpfr_mul(klpi, klpi, tmp2, MPFR_RNDN);
                        }
                        mpfr_div(klpi, klpi, fact->data[p], MPFR_RNDN);


                        mpfr_mul(skpn[idx + n], r1k, af, MPFR_RNDN);
                        mpfr_mul(skpn[idx + n], skpn[idx + n], klpi,
                                 MPFR_RNDN);
                        /* When compiling without the OpenMP parallelization,
                         * we could avoid the "mpfr_pow_ui" function by
                         * computing "tpow" recursively.  The recursion should
                         * be more efficient but does not allow parallelization
                         * over the "n" loop. */
                        mpfr_pow_ui(tpow, t, n + k + u + 1, MPFR_RNDN);
                        mpfr_mul(skpn[idx + n], skpn[idx + n], tpow,
                                 MPFR_RNDN);
                        mpfr_mul_ui(skpn[idx + n], skpn[idx + n], 2,
                                    MPFR_RNDN);
                        mpfr_div_ui(skpn[idx + n], skpn[idx + n], 2 * n + 1,
                                    MPFR_RNDN);
                    }
                    }


                    mpfr_clears(af, klpi, tpow, tmp1, tmp2, (mpfr_ptr)NULL);
                    mpfr_free_cache();
                }
            }
        }
    }


    mpfr_clears(r1k, rrefpow, t, (mpfr_ptr)NULL);
    CHARM(mpfr_ndarray_free)(fact);


    /* Some low-degree coefficients may not be defined, depending on "type", so
     * let's set them to "0.0" */
    for (unsigned k = kmin; k <= kmax; k++)
    {
        for (unsigned p = 1; p <= pmax; p++)
        {
            for (unsigned i = 0; i <= imax; i++)
            {
                idx = CHARM(misc_idx_4d)(k - kmin, p - 1, i, 0, pmax, imax + 1,
                                         nmax + 1);
                for (unsigned n = 0; n < CHARM_MIN(v, nmax + 1); n++)
                    mpfr_set_ui(skpn[idx + n], 0, MPFR_RNDN);
            }
        }
    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
EXIT:
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
    /* --------------------------------------------------------------------- */
}

