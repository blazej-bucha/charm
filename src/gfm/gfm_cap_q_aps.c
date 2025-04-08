/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../mpfr/mpfr_ndarray.h"
#include "../mpfr/mpfr_ndarray_check.h"
#include "gfm_cap_q_aps.h"
/* ------------------------------------------------------------------------- */






/* Computes the "aps" coefficients from Eq. (42) of Bucha et al. (2019a) */
void CHARM(gfm_cap_q_aps)(unsigned p,
                          unsigned s,
                          const mpfr_ndarray *fact,
                          mpfr_t aps,
                          CHARM(err) *err)
{
    if (p < 3)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"p\" must be larger than \"2\".");
        return;
    }


    if (CHARM(mpfr_ndarray_check)(fact, 1, (size_t)p))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"fact\" mpfr_ndarray.");
        return;
    }


    mpfr_mul(aps, fact->data[p - 1], fact->data[p - 3], MPFR_RNDN);
    mpfr_div(aps, aps, fact->data[p - s], MPFR_RNDN);
    mpfr_div(aps, aps, fact->data[p - s - 2], MPFR_RNDN);
    mpfr_div(aps, aps, fact->data[s - 1], MPFR_RNDN);
    if ((p - 1) % 2)
        mpfr_neg(aps, aps, MPFR_RNDN);


    mpfr_free_cache();


    return;
}

