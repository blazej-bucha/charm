/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../misc/misc_is_nearly_equal.h"
/* ------------------------------------------------------------------------- */






/* An internal function to check whether or not FFT can be applied for cell
 * grid-wise synthesis. */
_Bool CHARM(shs_grd_cell_fft_check)(const CHARM(cell) *grd,
                                    unsigned long nmax)
{
    if (grd->nlon > 1)
        /* Let's check whether FFT can be employed.  Six conditions must be
         * satisfied to allow FFT:
         *
         * i) the number of grid longitudes must be large enough when
         *    compared with "nmax",
         *
         * ii) the longitudinal step must be constant,
         *
         * iii) "grd->lonmin[0]" must be zero,
         *
         * iv) "grd->lonmax[grd->nlon - 1]" must be equal to "2.0 * PI",
         *
         * v) there must be at least two cells in the longitudinal
         *    direction, and
         *
         * vi) the maximum longitude of the "j"th cell must be equal to the
         *     minimum longitude of the "j + 1"th cell.
         *
         * The second condition has already been checked, so let's do the
         * remaining checks.  Due to the previous check, the sixth
         * condition needs to be checked only for the first to longitude
         * cells.*/
        if ((grd->nlon - 1) / 2 >= nmax &&
            CHARM(misc_is_nearly_equal)(grd->lonmin[0], PREC(0.0),
                                        CHARM(glob_threshold)) &&
            CHARM(misc_is_nearly_equal)(grd->lonmax[grd->nlon - 1],
                                        PREC(2.0) * PI,
                                        CHARM(glob_threshold)) &&
            CHARM(misc_is_nearly_equal)(grd->lonmin[1], grd->lonmax[0],
                                        CHARM(glob_threshold)))
            /* Great, FFT can be applied for this grid! */
            return 1;
        else
            /* Oh no, FFT cannot be applied for this grid. */
            return 0;
    else
        /* There is only one cell in the longitudinal direction, so no FFT
         * for simplicity.  In case, the speed is not at all crucial. */
        return 0;
}
