/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../prec.h"
#include "gfm_cap_q_norm_dist.h"
/* ------------------------------------------------------------------------- */






void CHARM(gfm_cap_q_norm_dist)(const mpfr_t t,
                                const mpfr_t t2,
                                const mpfr_t u,
                                mpfr_t distance)
{
    mpfr_mul(distance, t, u, MPFR_RNDN);
    mpfr_mul_si(distance, distance, -2, MPFR_RNDN);
    mpfr_add_ui(distance, distance, 1, MPFR_RNDN);
    mpfr_add(distance, distance, t2, MPFR_RNDN);
    mpfr_sqrt(distance, distance, MPFR_RNDN);


    mpfr_free_cache();


    return;
}

