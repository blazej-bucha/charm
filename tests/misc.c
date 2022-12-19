/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/prec.h"
#include "../src/simd/simd.h"
#include "abs_r.h"
/* ------------------------------------------------------------------------- */






/* Tests for various features of CHarm */
int misc(void)
{
    int errnum = 0;


    /* "ABS_R" macro */
    errnum += abs_r();


    return errnum;
}

