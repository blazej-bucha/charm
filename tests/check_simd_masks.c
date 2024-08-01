/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../src/prec.h"
#include "../src/simd/simd.h"
#include "check_simd_masks.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants */
/* ------------------------------------------------------------------------- */
#undef WARN_ALL_TRUE
#define WARN_ALL_TRUE "\n\n        WARNING: All-true mask test didn't " \
                      "pass!\n"


#undef WARN_ALL_FALSE
#define WARN_ALL_FALSE "\n\n        WARNING: All-false mask test " \
                       "didn't pass!\n"


#undef WARN_MIXED_TRUE_MIXED_FALSE
#define WARN_MIXED_TRUE_MIXED_FALSE "\n\n        WARNING: Mixed-true, " \
                                 "mixed-false mask test didn't pass!\n"
/* ------------------------------------------------------------------------- */






/* Check the "MASK_TRUE_ALL" macro from "../src/simd/simd.h". */
long int check_simd_mask_true_all(void)
{
    long int e = 0;


#ifdef SIMD
    /* SIMD instructions enabled */


    const RI_SIMD zero = SET_ZERO_RI;
    if (!MASK_TRUE_ALL(CAST_RI2R(EQ_RI(zero, zero))))  /* "mask" is all true */
    {
        printf(WARN_ALL_TRUE);
        e += 1;
    }


    if (MASK_TRUE_ALL(CAST_RI2R(GT_RI(zero, zero))))  /* "mask" is all false */
    {
        printf(WARN_ALL_FALSE);
        e += 1;
    }


    /* Create several mask such that "mask = [1, 0, 0, ..., 0, 0]", "mask = [1,
     * 1, 0, ..., 0, 0], ..., "mask = [1, 1, 1, ..., 1, 0]", where the length
     * of the arrays is given by "SIMD_SIZE". */
    INT x[SIMD_SIZE];
    const RI_SIMD mone = SET1_RI(-1);
    for (int j = 0; j < SIMD_SIZE - 1; j++)
    {
        INT k;
        for (k = 0; k <= j; k++)
            x[k] = 0;

        for (; k < SIMD_SIZE; k++)
            x[k] = -k;

        if (MASK_TRUE_ALL(CAST_RI2R(GT_RI(LOAD_RI(x), mone))))
        {
            printf(WARN_MIXED_TRUE_MIXED_FALSE);
            e += 1;
        }
    }
#else
    /* No SIMD */


    if (MASK_TRUE_ALL(0))
    {
        printf(WARN_ALL_TRUE);
        e += 1;
    }


    if (!MASK_TRUE_ALL(1))
    {
        printf(WARN_ALL_FALSE);
        e += 1;
    }
#endif


    return e;
}






/* Check the "MASK_TRUE_ANY" macro from "../src/simd/simd.h". */
long int check_simd_mask_true_any(void)
{
    long int e = 0;


#ifdef SIMD
    /* SIMD instructions enabled */


    const REAL_SIMD zero = SET_ZERO_R;
    if (!MASK_TRUE_ANY(EQ_R(zero, zero)))  /* "mask" is all true */
    {
        printf(WARN_ALL_TRUE);
        e += 1;
    }


    const REAL_SIMD one = SET1_R(PREC(1.0));
    if (MASK_TRUE_ANY(GE_R(zero, one)))  /* "mask" is all false */
    {
        printf(WARN_ALL_FALSE);
        e += 1;
    }


    /* Create several mask such that "mask = [1, 0, 0, ..., 0, 0]", "mask = [1,
     * 1, 0, ..., 0, 0], ..., "mask = [1, 1, 1, ..., 1, 0]", where the length
     * of the arrays is given by "SIMD_SIZE". */
    REAL x[SIMD_SIZE];
    for (int j = 0; j < SIMD_SIZE - 1; j++)
    {
        int k;
        for (k = 0; k <= j; k++)
            x[k] = 0;

        for (; k < SIMD_SIZE; k++)
            x[k] = -k;

        if (!MASK_TRUE_ANY(GE_R(LOAD_R(x), zero)))
        {
            printf(WARN_MIXED_TRUE_MIXED_FALSE);
            e += 1;
        }
    }
#else
    /* No SIMD */


    if (MASK_TRUE_ANY(0))
    {
        printf(WARN_ALL_TRUE);
        e += 1;
    }


    if (!MASK_TRUE_ANY(1))
    {
        printf(WARN_ALL_FALSE);
        e += 1;
    }
#endif


    return e;
}
