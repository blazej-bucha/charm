/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "gfm_check_p.h"
/* ------------------------------------------------------------------------- */






/* Checks the value of integer topography power. */
void CHARM(gfm_check_p)(unsigned p, CHARM(err) *err)
{
    if (p < 1)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Topography power must be positive.");


    return;
}
