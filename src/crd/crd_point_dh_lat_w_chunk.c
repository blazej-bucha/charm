/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "../misc/misc_is_nearly_equal.h"
#include "crd_point_dh_lat_w_chunk.h"
/* ------------------------------------------------------------------------- */






/* Computes latitudes and integration weights of Driscoll--Healy grids for
 * a given ``nmax`` and saves them to ``dh``. */
void CHARM(crd_point_dh_lat_w_chunk)(CHARM(point) *dh,
                                     unsigned long nmax,
                                     size_t local_nlat,
                                     size_t local_0_start,
                                     size_t local_nlat_north)
{
    unsigned long L    = nmax + 1;
    REAL L_fp          = (REAL)L;
    REAL c             = PREC(2.0) / L_fp;
    REAL M_PI_2L       = PI / (PREC(2.0) * L_fp);
    REAL thold         = CHARM(glob_threshold);
    unsigned long imin = local_0_start;
    unsigned long imax = local_0_start + local_nlat_north;


#if HAVE_OPENMP
#pragma omp parallel default(none) shared(nmax, dh, c, L, M_PI_2L) \
shared(imin, imax, local_nlat, local_0_start, thold)
#endif
    {
    REAL sclti, sclti0, sclti1, cclti2, sclti2, w_tmp, clti;
    unsigned long north, south;


    unsigned long i;
#if HAVE_OPENMP
#pragma omp for
#endif
    for (i = imin; i < imax; i++)
    {
        /* The "i"th co-latitude */
        clti = M_PI_2L * (REAL)(i);


        north = i - imin;
        south = local_nlat - 1 - north;
        if (local_0_start == 0)
            south += 1;  /* Because of the missing south pole */


        /* Latitudes */
        dh->lat[north] = PI_2 - clti;


        /* Integration weights.  The direct computation of the Driscoll--Healy
         * weights is painfully slow for high degrees due to the multiple
         * evaluation of "sin((2 * k + 1) * clti)".  To speed things up, we
         * apply the Chebyshev recurrence to evaluate these terms instead of
         * their direct evaluation.  This makes the code a bit longer and less
         * intuitive, but it's worth it. */
        /* ----------------------------------------------------------------- */
        sclti  = SIN(clti);


        if (nmax == 0)
            dh->w[north] = c * sclti;
        else
        {
            sclti0 = PREC(0.0);
            sclti1 = sclti;
            cclti2 = PREC(2.0) * COS(clti);


            w_tmp = sclti1;
            for (unsigned long k = 1; k < L; k++)
            {
                sclti2 = cclti2 * sclti1 - sclti0;


                sclti0 = sclti1;
                sclti1 = sclti2;
                sclti2 = cclti2 * sclti1 - sclti0;


                w_tmp += sclti2 / (REAL)(2 * k + 1);


                sclti0 = sclti1;
                sclti1 = sclti2;
            }


            dh->w[north] = c * sclti * w_tmp;
        }
        /* ----------------------------------------------------------------- */


        /* Apply the symmetry property. */
        /* ----------------------------------------------------------------- */
        /* Driscoll--Healy grids do not have the south pole, so do not apply
         * the symmetry property at the north pole */
        if ((local_0_start == 0) && (i == 0))
            continue;


        /* Do not apply the symmetry property at the equator.  This is so as
         * not to change the sign of the near-zero latitude. */
        if (CHARM(misc_is_nearly_equal)(dh->lat[north], PREC(0.0), thold))
            continue;


        dh->w[south]   = dh->w[north];
        dh->lat[south] = -dh->lat[north];
        /* ----------------------------------------------------------------- */
    }
    }


    return;
}
