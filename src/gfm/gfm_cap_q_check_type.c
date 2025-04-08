/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "gfm_cap_q_check_type.h"
/* ------------------------------------------------------------------------- */






/* Checks for the supported type of the trunction coefficients */
void CHARM(gfm_cap_q_check_type)(unsigned type,
                                 CHARM(err) *err)
{
    if ((type != CHARM_GFM_Q00) &&
        (type != CHARM_GFM_Q10) &&
        (type != CHARM_GFM_Q11) &&
        (type != CHARM_GFM_Q20) &&
        (type != CHARM_GFM_Q21) &&
        (type != CHARM_GFM_Q22))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported type of truncation coefficients.");


    return;
}
