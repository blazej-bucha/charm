/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/simd/simd.h"
#include "../src/misc/misc_is_nearly_equal.h"
#include "check_simd_blend_r.h"
#include "cmp_arrays.h"
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






/* ------------------------------------------------------------------------- */
#if CHARM_FLOAT
#   define TRUE_BITS 0xFFFFFFFF
#else
#   define TRUE_BITS 0xFFFFFFFFFFFFFFFFLL
#endif
union ir { INT i; REAL r;};
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


    const REAL TRUE_VAL  = PREC(9999.0);
    const REAL FALSE_VAL = PREC(-9999.0);


    const REAL_SIMD TRUE  = SET1_R(TRUE_VAL);
    const REAL_SIMD FALSE = SET1_R(FALSE_VAL);


    const REAL_SIMD ZERO = SET_ZERO_RI;
    const REAL_SIMD ONE  = SET1_RI(1);


    REAL_SIMD xv = BLEND_R(FALSE, TRUE, EQ_R(ZERO, ZERO));  /* "mask" is all
                                                             * true */
    REAL x[SIMD_SIZE];
    STORE_R(&x[0], xv);


    REAL xref[SIMD_SIZE];
    for (size_t i = 0; i < SIMD_SIZE; i++)
        xref[i] = TRUE_VAL;


    if (cmp_arrays(x, xref, SIMD_SIZE, CHARM(glob_threshold)))
    {
        printf(WARN_ALL_TRUE);
        e += 1;
    }


    xv = BLEND_R(FALSE, TRUE, EQ_R(ONE, ZERO));  /* "mask" is all false */
    STORE_R(&x[0], xv);


    for (size_t i = 0; i < SIMD_SIZE; i++)
        xref[i] = FALSE_VAL;


    if (cmp_arrays(x, xref, SIMD_SIZE, CHARM(glob_threshold)))
    {
        printf(WARN_ALL_FALSE);
        e += 1;
    }


    REAL_SIMD maskv;
    REAL mask[SIMD_SIZE];
    union ir true_bits;  true_bits.i  = TRUE_BITS;
    union ir false_bits; false_bits.i = 0;
    for (size_t j = 0; j < SIMD_SIZE - 1; j++)
    {
        size_t k;
        for (k = 0; k <= j; k++)
        {
            mask[k] = false_bits.r;
            xref[k] = FALSE_VAL;
        }


        for (; k < SIMD_SIZE; k++)
        {
            mask[k] = true_bits.r;
            xref[k] = TRUE_VAL;
        }


        maskv = LOAD_R(mask);
        xv    = BLEND_R(FALSE, TRUE, maskv);
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
