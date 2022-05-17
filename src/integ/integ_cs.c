/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
#include "integ_cs.h"
/* ------------------------------------------------------------------------- */






void CHARM(integ_cs)(REAL u0, REAL du, size_t nu, REAL a1, REAL a2, REAL *s)
{
    if (CHARM(misc_is_nearly_equal)(a2, ADDP(0.0), CHARM(glob_threshold)))
    /* "a2 == ADDP(0.0)" within the "threshold" */

        memset(s, 0, nu * sizeof(REAL));

    else if (CHARM(misc_is_nearly_equal)(FABS(a1), FABS(a2), 
                                         CHARM(glob_threshold)))
    /* "FABS(a1) == FABS(a2) != ADDP(0.0)" within the "threshold" */
    {
        REAL u1   = u0 +       du;
        REAL u2   = u0 + ADDP(2.0) * du;
        REAL k2   = ADDP(2.0) * a2;
        REAL k2_4 = ADDP(0.25) / a2;
        REAL d    = ADDP(2.0) * COS(k2 * du);


        REAL cc0 = COS(k2 * u1) - COS(k2 * u0);
        s[0] = -k2_4 * cc0;


        if (nu == 1)
            return;


        REAL cc1 = COS(k2 * u2) - COS(k2 * u1);
        s[1] = -k2_4 * cc1;


        if (nu == 2)
            return;


        REAL cc2;
        for (size_t j = 2; j < nu; j++)
        {
            cc2 = d * cc1 - cc0;


            s[j] = -k2_4 * cc2;


            cc0 = cc1;
            cc1 = cc2;
        }
    }
    else
    {
        REAL u1   = u0 +       du;
        REAL u2   = u0 + ADDP(2.0) * du;
        REAL k1   = a2 - a1;
        REAL k2   = a2 + a1;
        REAL k1_2 = ADDP(2.0) * k1;
        REAL k2_2 = ADDP(2.0) * k2;
        REAL d1   = ADDP(2.0) * COS(k1 * du);
        REAL d2   = ADDP(2.0) * COS(k2 * du);


        REAL cc1_0 = COS(k1 * u1) - COS(k1 * u0);
        REAL cc2_0 = COS(k2 * u1) - COS(k2 * u0);
        s[0] = -cc1_0 / k1_2 - cc2_0 / k2_2;


        if (nu == 1)
            return;


        REAL cc1_1 = COS(k1 * u2) - COS(k1 * u1);
        REAL cc2_1 = COS(k2 * u2) - COS(k2 * u1);
        s[1] = -cc1_1 / k1_2 - cc2_1 / k2_2;


        if (nu == 2)
            return;


        REAL cc1_2;
        REAL cc2_2;
        for (size_t j = 2; j < nu; j++)
        {
            cc1_2 = d1 * cc1_1 - cc1_0;
            cc2_2 = d2 * cc2_1 - cc2_0;


            s[j] = -cc1_2 / k1_2 - cc2_2 / k2_2;


            cc1_0 = cc1_1;
            cc1_1 = cc1_2;
            cc2_0 = cc2_1;
            cc2_1 = cc2_2;
        }
    }


    return;
}
