/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/simd/simd.h"
#include "../src/simd/calloc_aligned.h"
#include "../src/simd/free_aligned.h"
#include "../src/misc/misc_is_nearly_equal.h"
#include "check_simd_sum_r.h"
#include "cmp_vals.h"
/* ------------------------------------------------------------------------- */






/* Check the "SUM_R" macro from "../src/simd/simd.h". */
long int check_simd_sum_r(void)
{
    long int e = 0;


    /* Test array "x = [0, 1, 2, ..., SIMD_SIZE - 1]" */
    REAL *x = (REAL *)CHARM(calloc_aligned)(SIMD_MEMALIGN, SIMD_SIZE,
                                            sizeof(REAL));
    for (int i = 0; i < SIMD_SIZE; i++)
        x[i] = (REAL)i;


    /* Reference sum */
    REAL ref = (SIMD_SIZE * (SIMD_SIZE - PREC(1.0))) / PREC(2.0);


    e += cmp_vals_real(SUM_R(LOAD_R(x)), ref, CHARM(glob_threshold));


    CHARM(free_aligned)(x);


    return e;
}
