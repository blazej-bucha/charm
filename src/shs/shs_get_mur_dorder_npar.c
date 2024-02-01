/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../err/err_propagate.h"
#include "shs_check_grads.h"
#include "shs_point_gradn.h"
#include "shs_max_npar.h"
/* ------------------------------------------------------------------------- */






/* Internal function to get:
 *
 * * the multiplicative factor "GM / R^dorder" ("mur") that appears before the
 *   double summation in the synthesis, where "i" is the order of the radial
 *   derivative and "dorder" is explained below, and
 *
 * * the order of the potential derivative ("dorder"; "0" for potential, "1"
 *   for first-order derivative(s), "2" for second-order derivative(s)), and
 *
 * * the number of parameters/quantities to be synthesized ("npar"). */
void CHARM(shs_get_mur_dorder_npar)(const CHARM(shc) *shcs,
                                    int dr,
                                    int dlat,
                                    int dlon,
                                    REAL *mur,
                                    unsigned *dorder,
                                    size_t *npar,
                                    CHARM(err) *err)
{
    CHARM(shs_check_grads)(dr, dlat, dlon, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    if ((dr == GRAD_0) && (dlat == GRAD_0) && (dlon == GRAD_0))
    {
        *dorder    = 0;
        *npar      = 1;
    }
    else if ((dr == GRAD_1) && (dlat == GRAD_1) && (dlon == GRAD_1))
    {
        *dorder    = 1;
        *npar      = 3;
    }
    else if ((dr == GRAD_2) && (dlat == GRAD_2) && (dlon == GRAD_2))
    {
        *dorder    = 2;
        *npar      = 6;
    }
    else
    {
        *dorder    = dr + dlat + dlon;
        *npar      = 1;
    }


    REAL rpow = shcs->r;
    for (int i = 1; i <= *dorder; i++)
        rpow *= shcs->r;
    *mur = shcs->mu / rpow;


    return;
}
