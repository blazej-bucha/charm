/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(leg_func_r_ri)(unsigned long nmax, REAL *r, REAL *ri)
{
    REAL w;


    r[0] = ri[0] = PREC(0.0);
    for (unsigned long m = 1; m <= (2 * nmax + 3); m++)
    {
        w     = SQRT((REAL)m);
        r[m]  = w;
        ri[m] = PREC(1.0) / w;
    }


    return;
}
