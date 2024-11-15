/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "shc_check_distribution.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(shc_check_distribution)(const CHARM(shc) *shcs,
                                    CHARM(err) *err)
{
#if HAVE_MPI
    if (shcs->distributed)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "This function requires non-distributed "
                       "\"charm" CHARM_SUFFIX "_shc\" structure.");
#endif


    return shcs->distributed;
}
