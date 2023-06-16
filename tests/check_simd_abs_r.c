/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../src/prec.h"
#include "../src/simd/simd.h"
/* ------------------------------------------------------------------------- */






/* Check the "ABS_R" macro from "../src/simd/simd.h". */
long int check_simd_abs_r(void)
{
    long int e = 0;


    /* Tested value */
    REAL_SIMD x = SET1_R(PREC(-1.5));


    /* Reference value to be obtained after applying the "ABS_R" macro on
     * "x". */
    REAL_SIMD x_ref = SET1_R(PREC(1.5));


    ABS_R_INIT;
    if (MOVEMASK(EQ_R(ABS_R(x), x_ref)) != SIMD_TRUE)
    {
        printf("        WARNING: The \"ABS_R\" macro does not work "
               "correctly with negative floating point numbers!\n");
        e += 1;
    }



    x = SET1_R(PREC(1.5));
    if (MOVEMASK(EQ_R(ABS_R(x), x_ref)) != SIMD_TRUE)
    {
        printf("        WARNING: The \"ABS_R\" macro does not work "
               "correctly with positive floating point numbers!\n");
        e += 1;
    }


    return e;
}
