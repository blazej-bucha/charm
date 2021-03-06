/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(xnum_xnorm)(REAL *x, int *ix)
{
    REAL w = FABS(*x);
    if (w >= BIGS)
    {
        (*x)  *= BIGI;
        (*ix) += 1;
    }
    else if (w < BIGSI)
    {
        (*x)  *= BIG;
        (*ix) -= 1;
    }


    return;
}
