/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






/* An internal function to check cell boundaries. */
void CHARM(crd_check_cells)(const CHARM(crd) *cell, CHARM(err) *err)
{
    for (size_t i = 0; i < cell->nlat; i++)
    {
        if (cell->lat[2 * i] <= cell->lat[2 * i + 1])
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "Minimum cell latitudes must be smaller than "
                           "maximum cell latitudes.");
            return;
        }
    }


    for (size_t j = 0; j < cell->nlon; j++)
    {
        if (cell->lon[2 * j] >= cell->lon[2 * j + 1])
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "Minimum cell longitudes must be smaller than "
                           "maximum cell longitudes.");
            return;
        }
    }


    return;
}

