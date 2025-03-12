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
#include "gfm_cap_q_norm_dist.h"
#include "gfm_cap_q_ddist.h"
/* ------------------------------------------------------------------------- */






/* Derivatives of the reciprocal value of the distance "l" with respect to "r"
 * (multiple times) and with respect to "psi" (only once) (Eq. 101 of Bucha et
 * al., 2019) */
void CHARM(gfm_cap_q_ddist)(mpfr_ndarray *dl,
                            const mpfr_t r,
                            const mpfr_t rref,
                            const mpfr_t psi,
                            const mpfr_ndarray *fact,
                            const mpfr_ndarray *double_fact,
                            unsigned dmax_r,
                            unsigned dmax_psi,
                            mpfr_prec_t NBITS,
                            CHARM(err) *err)
{
    /* Check inputs */
    /* --------------------------------------------------------------------- */
    /* The order of the derivative with respect to "psi" cannot be larger than
     * "1" */
    if (dmax_psi > 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"dmax_psi\" cannot be larger than \"1\".");
        return;
    }


    if (CHARM(mpfr_ndarray_check)(fact, 1, dmax_r + 1))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"fact\" mpfr_ndarray.");
        return;
    }


    if (CHARM(mpfr_ndarray_check)(double_fact, 1, 2 * dmax_r))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"double_fact\" "
                       "mpfr_ndarray.");
        return;
    }


    if (CHARM(mpfr_ndarray_check)(dl, 1, dmax_r + 1))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"dl\" mpfr_ndarray.");
        return;
    }
    /* --------------------------------------------------------------------- */


    /* Compute the distance "l" */
    /* --------------------------------------------------------------------- */
    mpfr_t l, t, t2, cospsi, sinpsi, tmp1, tmp2, tmp3;
    mpfr_inits2(NBITS, l, t, t2, cospsi, sinpsi, tmp1, tmp2, tmp3,
                (mpfr_ptr)NULL);
    mpfr_cos(cospsi, psi, MPFR_RNDN);
    mpfr_sin(sinpsi, psi, MPFR_RNDN);


    mpfr_div(t, rref, r, MPFR_RNDN);
    mpfr_mul(t2, t, t, MPFR_RNDN);
    CHARM(gfm_cap_q_norm_dist)(t, t2, cospsi, l);
    mpfr_mul(l, l, r, MPFR_RNDN);
    /* --------------------------------------------------------------------- */


    /* Zero-order derivative of "1 / l" */
    /* --------------------------------------------------------------------- */
    mpfr_ui_div(dl->data[0], 1, l, MPFR_RNDN);
    if (dmax_psi == 1)
    {
        mpfr_mul(dl->data[0], dl->data[0], r, MPFR_RNDN);
        mpfr_mul(dl->data[0], dl->data[0], rref, MPFR_RNDN);
        mpfr_mul(dl->data[0], dl->data[0], sinpsi, MPFR_RNDN);
        mpfr_div(dl->data[0], dl->data[0], l, MPFR_RNDN);
        mpfr_div(dl->data[0], dl->data[0], l, MPFR_RNDN);
        mpfr_neg(dl->data[0], dl->data[0], MPFR_RNDN);
    }


    if (dmax_r == 0)
        goto EXIT;
    /* --------------------------------------------------------------------- */


    /* Higher order derivatives */
    /* --------------------------------------------------------------------- */
    int sign;


    for (unsigned k = 1; k <= dmax_r; k++)
    {
        mpfr_set_ui(dl->data[k], 0, MPFR_RNDN);


        for (unsigned s = 0; s <= k; s++)
        {
            if (!((k + s) % 2))
            {
                mpfr_pow_ui(tmp2, l, k + s + 1, MPFR_RNDN);
                mpfr_mul(tmp3, rref, cospsi, MPFR_RNDN);
                mpfr_sub(tmp3, r, tmp3, MPFR_RNDN);
                if (dmax_psi == 0)
                    mpfr_pow_si(tmp1, tmp3, s, MPFR_RNDN);
                else if (dmax_psi == 1)
                {
                    mpfr_pow_si(tmp1, tmp3, (long)s - 1, MPFR_RNDN);
                    mpfr_mul(tmp1, tmp1, rref, MPFR_RNDN);
                    mpfr_mul(tmp1, tmp1, sinpsi, MPFR_RNDN);
                    mpfr_mul(tmp3, tmp3, r, MPFR_RNDN);
                    mpfr_mul_ui(tmp3, tmp3, k + s + 1, MPFR_RNDN);
                    mpfr_div(tmp3, tmp3, l, MPFR_RNDN);
                    mpfr_div(tmp3, tmp3, l, MPFR_RNDN);
                    mpfr_ui_sub(tmp3, s, tmp3, MPFR_RNDN);
                    mpfr_mul(tmp1, tmp1, tmp3, MPFR_RNDN);
                }
                mpfr_div(tmp1, tmp1, tmp2, MPFR_RNDN);

                mpfr_mul(tmp1, tmp1, double_fact->data[k - s + 1], MPFR_RNDN);
                mpfr_mul(tmp1, tmp1, double_fact->data[k + s - 1], MPFR_RNDN);
                mpfr_div(tmp1, tmp1, fact->data[k - s + 1], MPFR_RNDN);
                mpfr_mul(tmp1, tmp1, fact->data[k], MPFR_RNDN);
                mpfr_div(tmp1, tmp1, fact->data[s], MPFR_RNDN);

                sign = (((k + s) / 2) % 2) ? -1 : 1;
                mpfr_mul_si(tmp1, tmp1, sign, MPFR_RNDN);

                mpfr_add(dl->data[k], dl->data[k], tmp1, MPFR_RNDN);
            }
        }

    }
    /* --------------------------------------------------------------------- */


    /* --------------------------------------------------------------------- */
EXIT:
    mpfr_clears(l, t, t2, cospsi, sinpsi, tmp1, tmp2, tmp3, (mpfr_ptr)NULL);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
    /* --------------------------------------------------------------------- */
}

