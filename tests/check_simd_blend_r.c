/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/simd/simd.h"
#include "cmp_arrays.h"
#include "check_simd_blend_r.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants */
/* ------------------------------------------------------------------------- */
#undef WARN_ALL_TRUE
#define WARN_ALL_TRUE "\n\n        WARNING: All-true blend test didn't " \
                      "pass!\n"


#undef WARN_ALL_FALSE
#define WARN_ALL_FALSE "\n\n        WARNING: All-false blend test " \
                       "didn't pass!\n"


#undef WARN_MIXED_TRUE_MIXED_FALSE
#define WARN_MIXED_TRUE_MIXED_FALSE "\n\n        WARNING: Mixed-true, " \
                                 "mixed-false blend test didn't pass!\n"
/* ------------------------------------------------------------------------- */






/* Check the "BLEND_R" macro from "../src/simd/simd.h".  If compiled without
 * the SIMD instructions, the function returns zero and is not called from the
 * main check program.  To make the program independent, it intentionally does
 * not use "MASK_TRUE_ALL" and "MASK_TRUE_ANY" macros. */
long int check_simd_blend_r(void)
{
#ifdef SIMD
    /* SIMD enabled */


    long int e = 0;


    const int TRUE_VAL_int  =  9999;
    const int FALSE_VAL_int = -9999;


    const REAL TRUE_VAL_REAL  = (REAL)TRUE_VAL_int;
    const REAL FALSE_VAL_REAL = (REAL)FALSE_VAL_int;


    const REAL_SIMD TRUE  = SET1_R(TRUE_VAL_REAL);
    const REAL_SIMD FALSE = SET1_R(FALSE_VAL_REAL);


    const REAL_SIMD ZERO = SET_ZERO_R;
    const REAL_SIMD ONE  = SET1_R(PREC(1.0));


    REAL_SIMD xv = BLEND_R(FALSE, TRUE, EQ_R(ZERO, ZERO));  /* "mask" is all
                                                             * true */
    REAL x[SIMD_SIZE];
    STORE_R(&x[0], xv);


    REAL xref[SIMD_SIZE];
    for (size_t i = 0; i < SIMD_SIZE; i++)
        xref[i] = TRUE_VAL_REAL;


    if (cmp_arrays(x, xref, SIMD_SIZE, CHARM(glob_threshold)))
    {
        printf(WARN_ALL_TRUE);
        e += 1;
    }


    xv = BLEND_R(FALSE, TRUE, EQ_R(ONE, ZERO));  /* "mask" is all false */
    STORE_R(&x[0], xv);


    for (size_t i = 0; i < SIMD_SIZE; i++)
        xref[i] = FALSE_VAL_REAL;


    if (cmp_arrays(x, xref, SIMD_SIZE, CHARM(glob_threshold)))
    {
        printf(WARN_ALL_FALSE);
        e += 1;
    }


    for (size_t j = 0; j < SIMD_SIZE - 1; j++)
    {
        size_t k;
        for (k = 0; k <= j; k++)
            xref[k] = FALSE_VAL_REAL;


        for (; k < SIMD_SIZE; k++)
            xref[k] = TRUE_VAL_REAL;


        xv = BLEND_R(FALSE, TRUE, EQ_R(TRUE, LOAD_R(xref)));
        STORE_R(&x[0], xv);


        if (cmp_arrays(x, xref, SIMD_SIZE, CHARM(glob_threshold)))
        {
            printf(WARN_MIXED_TRUE_MIXED_FALSE);
            e += 1;
        }
    }


    return e;
#else
    /* SIMD disabled */
    return 0;
#endif
}
