/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../src/mpfr/mpfr_is_nearly_equal.h"
#include "mpfr_cmp_vals.h"
/* ------------------------------------------------------------------------- */






/* Compares two real "mpfr_t" values up to some threshold */
/* ------------------------------------------------------------------------- */
long int mpfr_cmp_vals(mpfr_t x,
                       mpfr_t xref,
                       mpfr_t eps)
{
    long int ret = 0;


    if (!CHARM(mpfr_is_nearly_equal)(x, xref, eps))
    {
        printf("\n");
        printf("            result:    ");
        mpfr_out_str(stdout, 10, 0, x, MPFR_RNDN);
        printf("\n");
        printf("            reference: ");
        mpfr_out_str(stdout, 10, 0, xref, MPFR_RNDN);
        printf("\n");
        printf("            WARNING: Possibly too large error!\n");


        ret = 1;
    }


    return ret;
}
/* ------------------------------------------------------------------------- */

