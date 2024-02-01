/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../misc/misc_is_nearly_equal.h"
#include "shs_grd_point_fft_check.h"
/* ------------------------------------------------------------------------- */






/* An internal function to check whether or not FFT can be applied for point
 * grid-wise synthesis. */
_Bool CHARM(shs_grd_point_fft_check)(const CHARM(point) *grd,
                                     REAL deltalon,
                                     unsigned long nmax)
{
    /* Let's check whether FFT can be employed.  Four conditions must be
     * satisfied to allow FFT: i) the number of grid longitudes must be large
     * enough when compared with "nmax", ii) the longitudinal step must be
     * constant, iii) "grd->lon[0]" must be zero, and iv) "grd->lon[grd->nlon
     * 1] + deltalon" must be equal to "2.0 * PI".
     *
     * The second condition has already been checked, so let's do the
     * remaining
     * checks. */
    if ((grd->nlon - 1) / 2 >= nmax &&
        CHARM(misc_is_nearly_equal)(grd->lon[0], PREC(0.0),
                                    CHARM(glob_threshold)) &&
        CHARM(misc_is_nearly_equal)(grd->lon[grd->nlon - 1] + deltalon,
                                    PREC(2.0) * PI,
                                    CHARM(glob_threshold)))
        /* Great, FFT can be applied for this grid! */
        return 1;
    else
        /* Oh no, FFT cannot be applied for this grid. */
        return 0;
}
