/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "gfm_cap_q_check_radius.h"
/* ------------------------------------------------------------------------- */






/* Checks the spherical radius in cap-modified SGFM */
void CHARM(gfm_cap_q_check_radius)(const mpfr_t r,
                                   mpfr_prec_t NBITS,
                                   CHARM(err) *err)
{
    mpfr_t zero;
    mpfr_init2(zero, NBITS);
    mpfr_set_ui(zero, 0, MPFR_RNDN);


    if (mpfr_lessequal_p(r, zero))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Spherical radius must be positive.");


    mpfr_clear(zero);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
}

