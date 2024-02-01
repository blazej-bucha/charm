/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "shs_check_single_derivative.h"
#include "shs_point_gradn.h"
/* ------------------------------------------------------------------------- */






/* Internal function to check the orders of the radial, latitudinal and
 * longitudinal derivatives of the gravitational potential, including the full
 * gradients. */
void CHARM(shs_check_grads)(int dr,
                            int dlat,
                            int dlon,
                            CHARM(err) *err)
{
    if ((dr == GRAD_0) && (dlat == GRAD_0) && (dlon == GRAD_0))
        goto EXIT;
    else if ((dr == GRAD_1) && (dlat == GRAD_1) && (dlon == GRAD_1))
        goto EXIT;
    else if ((dr == GRAD_2) && (dlat == GRAD_2) && (dlon == GRAD_2))
        goto EXIT;


    CHARM(shs_check_single_derivative)(dr, dlat, dlon, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }


EXIT:
    return;
}
