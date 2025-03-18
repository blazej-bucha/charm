/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../mpfr/mpfr_ndarray.h"
#include "../mpfr/mpfr_ndarray_free.h"
#include "../mpfr/mpfr_ndarray_malloc.h"
#include "../mpfr/mpfr_ndarray_check.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "../misc/misc_idx_4d.h"
#include "gfm_cap_q_ddist.h"
#include "gfm_cap_q_aps.h"
#include "gfm_cap_q_dkernel.h"
/* ------------------------------------------------------------------------- */






void CHARM(gfm_cap_q_dkernel)(mpfr_ndarray *dk,
                              const mpfr_t r,
                              const mpfr_t rref,
                              const mpfr_t psi,
                              const mpfr_ndarray *rwq,
                              const mpfr_ndarray *binomial,
                              unsigned pmax,
                              unsigned kmax,
                              unsigned imax,
                              unsigned psimax,
                              const mpfr_ndarray *fact,
                              const mpfr_ndarray *double_fact,
                              mpfr_prec_t NBITS,
                              CHARM(err) *err)
{
    /* Check inputs */
    /* --------------------------------------------------------------------- */
    if (CHARM(mpfr_ndarray_check)(dk, 3, (size_t)pmax, (size_t)(imax + 1),
                                  (size_t)(kmax + 1)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"dk\" mpfr_ndarray.");
        return;
    }


    size_t nbin = CHARM_MAX(kmax + 1, imax + 1);
    if (CHARM(mpfr_ndarray_check)(binomial, 2, (size_t)nbin, (size_t)nbin))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"binomial\" mpfr_ndarray.");
        return;
    }


    if (CHARM(mpfr_ndarray_check)(fact, 1, (size_t)(pmax + 1)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"fact\" mpfr_ndarray.");
        return;
    }


    unsigned dmax = pmax + kmax + 1;
    if (CHARM(mpfr_ndarray_check)(double_fact, 1, (size_t)(2 * dmax)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"double_fact\" "
                       "mpfr_ndarray.");
        return;
    }


    if (CHARM(mpfr_ndarray_check)(rwq, 2, (size_t)(kmax + 1), (size_t)pmax))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"rwq\" mpfr_ndarray.");
        return;
    }
    /* --------------------------------------------------------------------- */


    /* Inits */
    /* --------------------------------------------------------------------- */
    mpfr_t tmp;
    mpfr_init2(tmp, NBITS);
    mpfr_ndarray *dl  = NULL;
    mpfr_ndarray *dk1 = NULL;
    /* --------------------------------------------------------------------- */


    /* Distance "l" and its derivatives "dl" (Eq. 68 of Bucha et al., 2019b) */
    /* --------------------------------------------------------------------- */
    /* Derivatives of "1 / l(r, psi)" with respect to "r" */
    dl = CHARM(mpfr_ndarray_malloc)(NBITS, 1, (size_t)(dmax + 1));
    if (dl == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    CHARM(gfm_cap_q_ddist)(dl, r, rref, psi, fact, double_fact, dmax, psimax,
                           NBITS, err);
    /* --------------------------------------------------------------------- */


    /* "i == 0" */
    /* ===================================================================== */
    /* K1 */
    /* --------------------------------------------------------------------- */
    dk1 = CHARM(mpfr_ndarray_malloc)(NBITS, 1, (size_t)(dmax + 1));
    if (dk1 == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EMEM,
                       CHARM_ERR_MALLOC_FAILURE);
        goto EXIT;
    }


    for (unsigned i = 0; i <= dmax; i++)
        mpfr_mul(dk1->data[i], rref, dl->data[i], MPFR_RNDN);
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
    for (size_t v = 0; v < dk->size; v++)
        mpfr_set_ui(dk->data[v], 0, MPFR_RNDN);
    /* --------------------------------------------------------------------- */


    /* p = 1 */
    /* --------------------------------------------------------------------- */
    size_t idx;
    for (unsigned k = 0; k <= kmax; k++)
    {
        idx = CHARM(misc_idx_4d)(0, 0, 0, k, pmax, imax + 1, kmax + 1);
        mpfr_set(dk->data[idx], dk1->data[k], MPFR_RNDN);
    }


    if (pmax == 1)
        goto IGT0;
    /* --------------------------------------------------------------------- */


    /* p = 2 */
    /* --------------------------------------------------------------------- */
    for (unsigned k = 0; k <= kmax; k++)
    {
        idx = CHARM(misc_idx_4d)(0, 1, 0, k, pmax, imax + 1, kmax + 1);


        mpfr_mul_si(dk->data[idx], dk1->data[k], -((long)k - 1), MPFR_RNDN);
        mpfr_mul(tmp, r, dk1->data[k + 1], MPFR_RNDN);
        mpfr_sub(dk->data[idx], dk->data[idx], tmp, MPFR_RNDN);
        mpfr_div_ui(dk->data[idx], dk->data[idx], 2, MPFR_RNDN);
    }


    if (pmax == 2)
        goto IGT0;
    /* --------------------------------------------------------------------- */


    /* p = 3 */
    /* --------------------------------------------------------------------- */
    mpfr_t aps, sum;
    mpfr_inits2(NBITS, aps, sum, (mpfr_ptr)NULL);


    for (unsigned k = 0; k <= kmax; k++)
    {
        for (unsigned p = 3; p <= pmax; p++)
        {
            idx = CHARM(misc_idx_4d)(0, p - 1, 0, k, pmax, imax + 1, kmax + 1);
            mpfr_set_ui(dk->data[idx], 0, MPFR_RNDN);


            for (unsigned s = 1; s <= p - 2; s++)
            {
                mpfr_set_ui(sum, 0, MPFR_RNDN);
                for (unsigned q = 0; q <= k; q++)
                {
                    mpfr_mul(tmp, binomial->data[k * binomial->shape[1] + q],
                             rwq->data[(k - q) * rwq->shape[1] + (p - s - 1)],
                             MPFR_RNDN);
                    mpfr_mul(tmp, tmp, dk1->data[p - s + q], MPFR_RNDN);
                    mpfr_add(sum, sum, tmp, MPFR_RNDN);
                }


                CHARM(gfm_cap_q_aps)(p, s, fact, aps, err);
                mpfr_mul(sum, sum, aps, MPFR_RNDN);
                mpfr_add(dk->data[idx], dk->data[idx], sum, MPFR_RNDN);
            }


            mpfr_div(dk->data[idx], dk->data[idx], fact->data[p], MPFR_RNDN);
        }
    }


    mpfr_clears(aps, sum, (mpfr_ptr)NULL);
    /* --------------------------------------------------------------------- */
    /* ===================================================================== */


    /* "i > 0" */
    /* ===================================================================== */
    size_t idxi0;
    unsigned smin;


IGT0:
    for (unsigned k = 0; k <= kmax; k++)
    {
        for (unsigned p = 1; p <= pmax; p++)
        {
            for (unsigned i = 1; i <= imax; i++)
            {
                idx = CHARM(misc_idx_4d)(0, p - 1, i, k, pmax, imax + 1,
                                         kmax + 1);


                /* "CHARM_MAX(1, p - i)", but in a special way, as "p - i" may
                 * be negative and "p" and "i" are unsigned ints. */
                smin = (p < i) ? 1: CHARM_MAX(1, p - i);
                for (unsigned s = smin; s <= p; s++)
                {
                    mpfr_set(tmp,
                             binomial->data[i * binomial->shape[1] + (p - s)],
                             MPFR_RNDN);
                    mpfr_mul_ui(tmp, tmp, s, MPFR_RNDN);
                    mpfr_div_ui(tmp, tmp, p, MPFR_RNDN);


                    idxi0 = CHARM(misc_idx_4d)(0, s - 1, 0, k, pmax, imax + 1,
                                               kmax + 1);
                    mpfr_mul(tmp, tmp, dk->data[idxi0], MPFR_RNDN);
                    mpfr_add(dk->data[idx], dk->data[idx], tmp, MPFR_RNDN);
                }
            }
        }
    }
    /* ===================================================================== */


EXIT:
    /* Clean up */
    /* --------------------------------------------------------------------- */
    mpfr_clear(tmp);
    CHARM(mpfr_ndarray_free)(dl);
    CHARM(mpfr_ndarray_free)(dk1);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;
    /* --------------------------------------------------------------------- */


    return;
}

