/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "mpfr_check_bits.h"
/* ------------------------------------------------------------------------- */






/* Internal function to check whether "NBITS" can be used to initialize
 * "mpfr_t". */
void CHARM(mpfr_check_bits)(mpfr_prec_t NBITS, CHARM(err) *err)
{
    if (NBITS > MPFR_PREC_MAX)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Too large value for the \"mpfr_prec_t\" data type.");
    else if (NBITS < MPFR_PREC_MIN)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Too small value for the \"mpfr_prec_t\" data type.");


    return;
}

