/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../mpfr/mpfr_ndarray.h"
#include "../mpfr/mpfr_ndarray_check.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "gfm_cap_q_dr.h"
/* ------------------------------------------------------------------------- */






/* Eq. 75 of Bucha et al. (2019) */
void CHARM(gfm_cap_q_dr)(mpfr_ndarray *rps,
                         const mpfr_ndarray *rpows,
                         unsigned pmax,
                         unsigned kmax,
                         mpfr_prec_t NBITS,
                         CHARM(err) *err)
{
    /* Check inputs */
    /* --------------------------------------------------------------------- */
    if (CHARM(mpfr_ndarray_check)(rps, 2, (size_t)(kmax + 1), (size_t)pmax))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"rps\" mpfr_ndarray.");
        return;
    }


    if (CHARM(mpfr_ndarray_check)(rpows, 1, (size_t)(pmax + 1)))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"rpows\" mpfr_ndarray.");
        return;
    }
    /* --------------------------------------------------------------------- */


    size_t idx;
    mpfr_t tmp;
    mpfr_init2(tmp, NBITS);


    for (unsigned k = 0; k <= kmax; k++)
    {
        for (unsigned p = 1; p <= pmax; p++)
        {
            idx = k * pmax + p - 1;


            if (k == 0)
                mpfr_set(rps->data[idx], rpows->data[p], MPFR_RNDN);
            else if (k > p)
                mpfr_set_ui(rps->data[idx], 0, MPFR_RNDN);
            else
            {
                mpfr_set_ui(rps->data[idx], 1, MPFR_RNDN);
                for (unsigned j = 1; j <= k; j++)
                    mpfr_mul_ui(rps->data[idx], rps->data[idx], p - j + 1,
                                MPFR_RNDN);
                mpfr_mul(rps->data[idx], rps->data[idx], rpows->data[p - k],
                         MPFR_RNDN);
            }
        }
    }


    mpfr_clear(tmp);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
}

