/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_propagate.h"
#include "../err/err_set.h"
#include "../misc/misc_is_nearly_equal.h"
#include "../misc/misc_arr_chck_lin_incr.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_check_grd_lons)(const CHARM(crd) *cell, REAL *dlon,
                                    CHARM(err) *err)
{
    if (cell->nlon > 1)
    {
        int err_tmp = CHARM(misc_arr_chck_lin_incr)(cell->lon, 2 * cell->nlon,
                                                    0, 2,
                                                    CHARM(glob_threshold2),
                                                    err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
        if (err_tmp != 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "\"cell->lon\" is not a linearly increasing array "
                           "of cells within the \"threshold2\".");
            return;
        }


        err_tmp = CHARM(misc_arr_chck_lin_incr)(cell->lon, 2 * cell->nlon,
                                                1, 2, CHARM(glob_threshold2),
                                                err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
        if (err_tmp != 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "\"cell->lon\" is not a linearly increasing array "
                           "of cells within the \"threshold2\".");
            return;
        }
    }


    /* At this point, we know that the steps between cells in "cell->lon" are
     * constant.  Here, we check if the steps between the minimum longitudes
     * and the maximum longitudes are equal and, if true, store this value in
     * a separate variable (will be necessary later for the PSLR algorithm) */
    if (cell->nlon > 1)
    {
        if (CHARM(misc_is_nearly_equal)(cell->lon[2] - cell->lon[0],
                                        cell->lon[3] - cell->lon[1],
                                        CHARM(glob_threshold2)) == 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The difference \"cell->lon[2] - cell->lon[0]\" "
                           "has to be equal to "
                           "\"cell->lon[3] - cell->lon[1]\".");
            return;
        }
        *dlon = cell->lon[2] - cell->lon[0];
    }
    else
        *dlon = cell->lon[1] - cell->lon[0];


    return;
}
