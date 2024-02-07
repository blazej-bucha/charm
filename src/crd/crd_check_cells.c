/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "crd_check_cells.h"
/* ------------------------------------------------------------------------- */






/* An internal function to check cell boundaries. */
void CHARM(crd_check_cells)(const CHARM(cell) *cell, CHARM(err) *err)
{
    for (size_t i = 0; i < cell->nlat; i++)
    {
        if (cell->latmax[i] <= cell->latmin[i])
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "Minimum cell latitudes must be smaller than "
                           "maximum cell latitudes.");
            return;
        }
    }


    for (size_t j = 0; j < cell->nlon; j++)
    {
        if (cell->lonmin[j] >= cell->lonmax[j])
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "Minimum cell longitudes must be smaller than "
                           "maximum cell longitudes.");
            return;
        }
    }


    return;
}

