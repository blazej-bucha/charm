/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "check_simd_abs_r.h"
#include "check_simd_neg_r.h"
#include "check_func.h"
#include "check_outcome.h"
#include "misc.h"
/* ------------------------------------------------------------------------- */






/* Tests for various features of CHarm */
long int misc(void)
{
    long int e    = 0;
    long int esum = 0;


    check_func("ABS_R");
    e = check_simd_abs_r();
    check_outcome(e);
    esum += e;


    check_func("NEG_R");
    e = check_simd_neg_r();
    check_outcome(e);
    esum += e;


    return esum;
}

