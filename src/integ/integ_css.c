/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "integ_ss.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(integ_css)(REAL u0, REAL du, REAL a1, REAL a2)
{
    REAL i1, i2;
    CHARM(integ_ss)(u0, du, 1, a1 + a2, ADDP(1.0), &i1);
    CHARM(integ_ss)(u0, du, 1, a1 - a2, ADDP(1.0), &i2);


    return (i1 - i2) / ADDP(2.0);
}
