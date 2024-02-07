/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include "xnum_x2f.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(xnum_x2f)(REAL x, int ix)
{
    if (ix == 0)
        return x;
    else if (ix < -1)
        return PREC(0.0);
    else if (ix < 0)
        return (x * BIGI);
    else
        return (x * BIG);
}
