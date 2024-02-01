/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <math.h>
#include "../prec.h"
#include "xnum_xnorm.h"
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
