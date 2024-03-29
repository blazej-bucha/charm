/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "misc_is_nearly_equal.h"
/* ------------------------------------------------------------------------- */






_Bool CHARM(misc_is_nearly_equal)(REAL a, REAL b, REAL eps)
{
    REAL a_abs    = FABS(a);
    REAL b_abs    = FABS(b);
    REAL diff_abs = FABS(a - b);
    REAL sum_abs  = a_abs + b_abs;


    if (a == b)
        return 1;
    else if (a == PREC(0.0) || b == PREC(0.0))
        return (diff_abs <= eps);
    else if (sum_abs == a_abs || sum_abs == b_abs)
        return (diff_abs <= eps);
    else
        return (diff_abs / CHARM_MAX(a_abs, b_abs)) <= eps;
}
