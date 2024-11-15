/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "err_check_distribution.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(err_check_distribution)(CHARM(err) *err)
{
#if HAVE_MPI
    if (err->distributed)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "This function requires non-distributed "
                       "\"charm" CHARM_SUFFIX "_err\" structure.");
#endif


    return err->distributed;
}
