/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/simd/simd.h"
#include "check_simd_abs_r.h"
#include "check_simd_neg_r.h"
#include "check_simd_sum_r.h"
#include "check_simd_masks.h"
#include "check_simd_blend_r.h"
#include "check_func.h"
#include "check_outcome.h"
#include "misc.h"
/* ------------------------------------------------------------------------- */






/* Tests for various features of CHarm */
long int misc(void)
{
    long int e    = 0;
    long int esum = 0;


    check_func("MASK_TRUE_ALL");
    e = check_simd_mask_true_all();
    check_outcome(e);
    esum += e;


    check_func("MASK_TRUE_ANY");
    e = check_simd_mask_true_any();
    check_outcome(e);
    esum += e;


    check_func("ABS_R");
    e = check_simd_abs_r();
    check_outcome(e);
    esum += e;


    check_func("NEG_R");
    e = check_simd_neg_r();
    check_outcome(e);
    esum += e;


    check_func("SUM_R");
    e = check_simd_sum_r();
    check_outcome(e);
    esum += e;


#ifdef SIMD
    check_func("BLEND_R");
    e = check_simd_blend_r();
    check_outcome(e);
    esum += e;
#endif


    return esum;
}

