/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include <stdio.h>
#include <mpfr.h>
#include "../err/err_set.h"
#include "mpfr_flush_unreleased_memory.h"
#include "mpfr_ndarray.h"
#include "mpfr_ndarray_check.h"
#include "mpfr_binomial.h"
/* ------------------------------------------------------------------------- */






/* Internal function to compute binomial coefficients for all integer
 * combinations up to "n" using MPFR.  The binomial coefficient "i" choose "k"
 * is obtained as "binom->data[i * binom->shape[1] + k]".  */
void CHARM(mpfr_binomial)(mpfr_ndarray *binom,
                          unsigned n,
                          mpfr_prec_t NBITS,
                          CHARM(err) *err)
{
    if (CHARM(mpfr_ndarray_check)(binom, 2, n + 1, n + 1))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong shape of the input \"binom\" mpfr_ndarray.");
        return;
    }


    /* Set all array elements of "binomial->data" to zero */
    for (size_t i = 0; i < binom->size; i++)
        mpfr_set_ui(binom->data[i], 0, MPFR_RNDN);


    mpfr_t tmp1;
    mpfr_init2(tmp1, NBITS);
    size_t idx;
    size_t offset = binom->shape[1];


    /* Compute the binomial coefficients */
    for (unsigned i = 0; i <= n; i++)
    {
        idx = i * offset;
        mpfr_set_ui(binom->data[idx], 1, MPFR_RNDN);


        for (unsigned j = 1; j <= i; j++)
        {
            mpfr_set_ui(tmp1, i - j + 1, MPFR_RNDN);
            mpfr_div_ui(tmp1, tmp1, j, MPFR_RNDN);
            mpfr_mul(binom->data[idx + j], binom->data[idx + j - 1], tmp1,
                     MPFR_RNDN);
        }
    }


    mpfr_clear(tmp1);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
}

