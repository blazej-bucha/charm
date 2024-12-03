/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "mpi_err_isdistributed.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(mpi_err_isdistributed)(CHARM(err) *err)
{
    if (!err->distributed)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The error structure \"charm" CHARM_SUFFIX "_err\" "
                       "must be initialized using \"charm" CHARM_SUFFIX
                       "_mpi_err_init\".");


    return err->distributed;
}
