/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "../prec.h"
#include "mpfr_flush_unreleased_memory.h"
#include "mpfr_is_nearly_equal.h"
/* ------------------------------------------------------------------------- */






_Bool CHARM(mpfr_is_nearly_equal)(mpfr_t a,
                                  mpfr_t b,
                                  mpfr_t eps)
{
    mpfr_t a_abs;
    mpfr_init2(a_abs, mpfr_get_prec(a));

    mpfr_t b_abs;
    mpfr_init2(b_abs, mpfr_get_prec(b));

    mpfr_t diff_abs;
    mpfr_prec_t abp = CHARM_MIN(mpfr_get_prec(a), mpfr_get_prec(b));
    mpfr_init2(diff_abs, abp);

    mpfr_t sum_abs;
    mpfr_init2(sum_abs, abp);

    mpfr_t tmp;
    mpfr_init2(tmp, abp);


    mpfr_abs(a_abs, a, MPFR_RNDN);
    mpfr_abs(b_abs, b, MPFR_RNDN);
    mpfr_sub(diff_abs, a, b, MPFR_RNDN);
    mpfr_abs(diff_abs, diff_abs, MPFR_RNDN);
    mpfr_add(sum_abs, a_abs, b_abs, MPFR_RNDN);


    _Bool ret = mpfr_cmp(a, b) == 0;
    if (ret)  /* "a == b" */
        goto EXIT;


    if (mpfr_zero_p(a) || mpfr_zero_p(b))  /* "(a == 0.0) || (b == 0.0)" */
    {
        ret = mpfr_cmp(diff_abs, eps) <= 0;
        goto EXIT;
    }


    if ((mpfr_cmp(sum_abs, a_abs) == 0) ||
        (mpfr_cmp(sum_abs, b_abs) == 0))  /* "(diff_abs == a_abs) ||
                                           *  (diff_abs == b_abs)" */
    {
        ret = mpfr_cmp(diff_abs, eps) <= 0;
        goto EXIT;
    }


    /* Else */
    mpfr_max(tmp, a_abs, b_abs, MPFR_RNDN);
    mpfr_div(tmp, diff_abs, tmp, MPFR_RNDN);
    ret = mpfr_cmp(tmp, eps) <= 0;


EXIT:
    mpfr_clears(a_abs, b_abs, diff_abs, sum_abs, tmp, (mpfr_ptr)NULL);
    FLUSH_UNRELEASED_MEMORY;
    return ret;
}
