/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "shs_check_single_derivative.h"
/* ------------------------------------------------------------------------- */






/* Internal function to check the orders of the radial, latitudinal and
 * longitudinal derivatives of the gravitational potential. */
void CHARM(shs_check_single_derivative)(int dr,
                                        int dlat,
                                        int dlon,
                                        CHARM(err) *err)
{
    char err_msg[CHARM_ERR_MAX_MSG];


    if (dr < 0)
    {
        snprintf(err_msg, CHARM_ERR_MAX_MSG,
                 "\"dr\" is \"%d\", but it must be non-negative.", dr);
        CHARM(err_set)(err, __FILE__, __LINE__,  __func__, CHARM_EFUNCARG,
                       err_msg);
        goto EXIT;
    }


    if (dlat < 0)
    {
        snprintf(err_msg, CHARM_ERR_MAX_MSG,
                 "\"dlat\" is \"%d\", but it must be non-negative.", dlat);
        CHARM(err_set)(err, __FILE__, __LINE__,  __func__, CHARM_EFUNCARG,
                       err_msg);
        goto EXIT;
    }


    if (dlon < 0)
    {
        snprintf(err_msg, CHARM_ERR_MAX_MSG,
                 "\"dlon\" is \"%d\", but it must be non-negative.", dlon);
        CHARM(err_set)(err, __FILE__, __LINE__,  __func__, CHARM_EFUNCARG,
                       err_msg);
        goto EXIT;
    }


    if ((dr + dlat + dlon) > SHS_MAX_DERIVATIVE)
    {
        snprintf(err_msg, CHARM_ERR_MAX_MSG,
                         "The sum \"dr + dlat + dlon\" is \"%d\", "
                         "but it cannot be larger than \"%d\".",
                         dr + dlat + dlon, SHS_MAX_DERIVATIVE);
        CHARM(err_set)(err, __FILE__, __LINE__,  __func__, CHARM_EFUNCARG,
                       err_msg);
        goto EXIT;
    }


EXIT:
    return;
}
