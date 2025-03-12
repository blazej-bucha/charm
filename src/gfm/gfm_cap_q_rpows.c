/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../mpfr/mpfr_ndarray.h"
#include "../mpfr/mpfr_ndarray_check.h"
#include "gfm_cap_q_rpows.h"
/* ------------------------------------------------------------------------- */






void CHARM(gfm_cap_q_rpows)(mpfr_ndarray *rpows,
                            const mpfr_t r,
                            CHARM(err) *err)
{
    if (rpows->ndim != 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong dimensions of the input \"rpows\" "
                       "mpfr_ndarray.");
        return;
    }


    mpfr_set_ui(rpows->data[0], 1, MPFR_RNDN);
    for (unsigned i = 1; i < rpows->shape[0]; i++)
        mpfr_mul(rpows->data[i], rpows->data[i - 1], r, MPFR_RNDN);


    mpfr_free_cache();


    return;
}

