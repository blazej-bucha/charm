/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "integ_cs.h"
#include "integ_sc.h"
/* ------------------------------------------------------------------------- */






void CHARM(integ_sc)(REAL u0, REAL du, size_t nu, REAL a1, REAL a2, REAL *s)
{
    /* Luckily, if we switch the "a1" and "a2" variables, we can take advantage
     * of the "CHARM(integ_cs)" function. */
    CHARM(integ_cs)(u0, du, nu, a2, a1, s);


    return;
}
