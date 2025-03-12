/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include <stdio.h>
#include <mpfr.h>
#include "mpfr_double_fact.h"
/* ------------------------------------------------------------------------- */






/* Internal function to compute the double factorial of "x" and store it in
 * "out". */
void CHARM(mpfr_double_fact)(unsigned x,
                             mpfr_t out)
{
    mpfr_set_ui(out, 1, MPFR_RNDN);
    for (unsigned i = x; i > 1; i -= 2)
        mpfr_mul_ui(out, out, i, MPFR_RNDN);


    mpfr_free_cache();


    return;
}
