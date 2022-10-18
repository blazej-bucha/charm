/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Computes latitudes and integration weights of Driscoll--Healy grids for
 * a given ``nmax`` and saves them to ``dh``. */
void CHARM(crd_point_dh_lats_weights)(CHARM(point) *dh, unsigned long nmax)
{
    unsigned long L  = nmax + 1;
    unsigned long L2 = 2 * L;
    REAL L_fp        = (REAL)L;
    REAL c           = PREC(2.0) / L_fp;
    REAL M_PI_2L     = PI / (PREC(2.0) * L_fp);


#if CHARM_PARALLEL
#pragma omp parallel default(none) shared(nmax, dh, c, L, L2, M_PI_2L)
#endif
    {
    REAL sclti, sclti0, sclti1, cclti2, sclti2, w_tmp, clti;


#if CHARM_PARALLEL
#pragma omp for
#endif
    for (unsigned long i = 0; i < L2; i++)
    {
        /* The "i"th co-latitude */
        clti = M_PI_2L * (REAL)(i);


        /* Latitudes */
        dh->lat[i] = PI_2 - clti;


        /* Integration weights.  The direct computation of the Driscoll--Healy
         * weights is painfully slow for high degrees due to the multiple
         * evaluation of "sin((2 * k + 1) * clti)".  To speed things up, we
         * apply the Chebyshev recurrence to evaluate these terms instead of
         * their direct evaluation.  This makes the code a bit longer and less
         * intuitive, but it's worth it. */
        /* ----------------------------------------------------------------- */
        sclti  = SIN(clti);


        if (nmax == 0)
            dh->w[i] = c * sclti;
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


            dh->w[i] = c * sclti * w_tmp;
        }
        /* ----------------------------------------------------------------- */
    }
    }


    return;
}
