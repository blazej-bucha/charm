/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "gfm_cap_q_check_psi.h"
/* ------------------------------------------------------------------------- */






/* Check the integration radius in cap-modified spectral gravity forward
 * modelling */
void CHARM(gfm_cap_q_check_psi)(const mpfr_t psi,
                                mpfr_prec_t NBITS,
                                CHARM(err) *err)
{
    mpfr_t pi, zero;
    mpfr_inits2(NBITS, pi, zero, (mpfr_ptr)NULL);
    mpfr_set_ui(zero, 0, MPFR_RNDN);
    mpfr_const_pi(pi, MPFR_RNDN);


    if (mpfr_less_p(psi, zero))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"psi\" cannot be smaller than zero.");
        goto EXIT;
    }


    if (mpfr_greater_p(psi, pi))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"psi\" cannot be larger than \"pi\".");
        goto EXIT;
    }


EXIT:
    mpfr_clears(pi, zero, (mpfr_ptr)NULL);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
}

