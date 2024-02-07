/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "integ_css.h"
#include "integ_scs.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(integ_scs)(REAL u0, REAL du, REAL a1, REAL a2)
{
    /* Luckily, if we switch the "a1" and "a2" variables, we can take advantage
     * of the "CHARM(integ_cs)" function. */
    return CHARM(integ_css)(u0, du, a2, a1);
}
