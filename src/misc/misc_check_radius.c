/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "misc_check_radius.h"
/* ------------------------------------------------------------------------- */






void CHARM(misc_check_radius)(REAL r,
                              CHARM(err) *err)
{
    if (r <= PREC(0.0))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Spherical radius must be positive.");


    return;
}
