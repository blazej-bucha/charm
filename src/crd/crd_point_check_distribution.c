/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "crd_point_check_distribution.h"
/* ------------------------------------------------------------------------- */







_Bool CHARM(crd_point_check_distribution)(const CHARM(point) *pnt,
                                          CHARM(err) *err)
{
#if HAVE_MPI
    if (pnt->distributed)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "This function requires non-distributed "
                       "\"charm" CHARM_SUFFIX "_point\" structure.");
#endif


    return pnt->distributed;
}
