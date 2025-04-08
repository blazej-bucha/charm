/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../prec.h"
#include "../mpfr/mpfr_ndarray.h"
#include "../mpfr/mpfr_ndarray_check.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "../err/err_set.h"
#include "gfm_cap_q_dnorm_dist.h"
/* ------------------------------------------------------------------------- */






/* Derivatives of the reciprocal value of the normalized distance (Eq. 61 of
 * Bucha et al., 2019) */
void CHARM(gfm_cap_q_dnorm_dist)(mpfr_ndarray *dg,
                                 const mpfr_t g,
                                 const mpfr_t t,
                                 const mpfr_t u,
                                 const mpfr_ndarray *fact,
                                 const mpfr_ndarray *double_fact,
                                 unsigned dmax,
                                 mpfr_prec_t NBITS,
                                 CHARM(err) *err)
{
    /* Check inputs */
    /* --------------------------------------------------------------------- */
    if (CHARM(mpfr_ndarray_check)(dg, 1, (size_t)(dmax + 1)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"dg\" mpfr_ndarray.");
        return;
    }


    mpfr_ui_div(dg->data[0], 1, g, MPFR_RNDN);
    if (dmax == 0)
        return;


    if (CHARM(mpfr_ndarray_check)(fact, 1, (size_t)(dmax + 1)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"fact\" mpfr_ndarray.");
        return;
    }


    if (CHARM(mpfr_ndarray_check)(double_fact, 1, (size_t)(2 * dmax + 1)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"double_fact\" "
                       "mpfr_ndarray.");
        return;
    }
    /* --------------------------------------------------------------------- */


    int sign;
    mpfr_t tmp1, tmp2;
    mpfr_inits2(NBITS, tmp1, tmp2, (mpfr_ptr)NULL);


    for (unsigned i = 1; i <= dmax; i++)
    {
        mpfr_set_ui(dg->data[i], 0, MPFR_RNDN);


        for (unsigned s = 0; s <= i; s++)
        {
            if (!((i + s) % 2))
            {
                mpfr_pow_ui(tmp2, g, i + s + 1, MPFR_RNDN);
                mpfr_sub(tmp1, t, u, MPFR_RNDN);
                mpfr_pow_ui(tmp1, tmp1, s, MPFR_RNDN);
                mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);

                mpfr_mul(tmp1, tmp1, double_fact->data[i - s + 1], MPFR_RNDN);
                mpfr_mul(tmp1, tmp1, double_fact->data[i + s - 1], MPFR_RNDN);
                mpfr_div(tmp1, tmp1, fact->data[i - s + 1], MPFR_RNDN);
                mpfr_mul(tmp1, tmp1, fact->data[i], MPFR_RNDN);
                mpfr_div(tmp1, tmp1, fact->data[s], MPFR_RNDN);

                sign = ((i + s) / 2) % 2 ? -1 : 1;
                mpfr_mul_si(tmp1, tmp1, sign, MPFR_RNDN);

                mpfr_add(dg->data[i], dg->data[i], tmp1, MPFR_RNDN);
            }
        }

    }


    mpfr_clears(tmp1, tmp2, (mpfr_ptr)NULL);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
}

