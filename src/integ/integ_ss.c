/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../prec.h"
#include "../misc/misc_is_nearly_equal.h"
/* ------------------------------------------------------------------------- */






void CHARM(integ_ss)(REAL u0, REAL du, size_t nu, REAL a1, REAL a2, REAL *s)
{
    REAL abs_a1 = FABS(a1);


    if (CHARM(misc_is_nearly_equal)(a1, PREC(0.0), CHARM(glob_threshold))
     || CHARM(misc_is_nearly_equal)(a2, PREC(0.0), CHARM(glob_threshold)))
    /* "(a1 == 0) || (a2 == 0)" within the "threshold" */

        memset(s, 0, nu * sizeof(REAL));

    else if (CHARM(misc_is_nearly_equal)(abs_a1, FABS(a2),
                                         CHARM(glob_threshold)))
    /* "FABS(a1) == FABS(a2) != PREC(0.0)" within the "threshold" */
    {
        /* We already know that "a1 != 0" and "a2 != 0", so that it is now
         * sufficient to check only whether "a1" and "a2" are positive or
         * negative. */
        REAL sign = ((a1 * a2) > 0) ? PREC(1.0) : PREC(-1.0);
        REAL u1   = u0 +       du;
        REAL u2   = u0 + PREC(2.0) * du;
        REAL du2  = PREC(0.5) * du;
        REAL k1   = PREC(2.0) * abs_a1;
        REAL a1_4 = PREC(0.25) / abs_a1;
        REAL d    = PREC(2.0) * COS(k1 * du);


        REAL ss0 = SIN(k1 * u1) - SIN(k1 * u0);
        s[0] = sign * (du2 - a1_4 * ss0);


        if (nu == 1)
            return;


        REAL ss1 = SIN(k1 * u2) - SIN(k1 * u1);
        s[1] = sign * (du2 - a1_4 * ss1);


        if (nu == 2)
            return;


        REAL ss2;
        for (size_t j = 2; j < nu; j++)
        {
            ss2 = d * ss1 - ss0;


            s[j] = sign * (du2 - a1_4 * ss2);


            ss0 = ss1;
            ss1 = ss2;
        }
    }
    else
    {
        REAL u1   = u0 +       du;
        REAL u2   = u0 + PREC(2.0) * du;
        REAL k1   = a2 - a1;
        REAL k2   = a2 + a1;
        REAL k1_2 = PREC(2.0) * k1;
        REAL k2_2 = PREC(2.0) * k2;
        REAL d1   = PREC(2.0) * COS(k1 * du);
        REAL d2   = PREC(2.0) * COS(k2 * du);


        REAL ss1_0 = SIN(k1 * u1) - SIN(k1 * u0);
        REAL ss2_0 = SIN(k2 * u1) - SIN(k2 * u0);
        s[0] = ss1_0 / k1_2 - ss2_0 / k2_2;


        if (nu == 1)
            return;


        REAL ss1_1 = SIN(k1 * u2) - SIN(k1 * u1);
        REAL ss2_1 = SIN(k2 * u2) - SIN(k2 * u1);
        s[1] = ss1_1 / k1_2 - ss2_1 / k2_2;


        if (nu == 2)
            return;


        REAL ss1_2;
        REAL ss2_2;
        for (size_t j = 2; j < nu; j++)
        {
            ss1_2 = d1 * ss1_1 - ss1_0;
            ss2_2 = d2 * ss2_1 - ss2_0;


            s[j] = ss1_2 / k1_2 - ss2_2 / k2_2;


            ss1_0 = ss1_1;
            ss1_1 = ss1_2;
            ss2_0 = ss2_1;
            ss2_1 = ss2_2;
        }
    }


    return;
}
