/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "mpfr_cmp_vals.h"
#include "mpfr_cmp_arrays.h"
/* ------------------------------------------------------------------------- */






/* Compares elements of two "mpfr_t" arrays up to some threshold */
/* ------------------------------------------------------------------------- */
long int mpfr_cmp_arrays(mpfr_t *arr1,
                         mpfr_t *arr2,
                         size_t n,
                         mpfr_t eps)
{
    long int ret = 0;


    for (size_t i = 0; i < n; i++)
        ret += mpfr_cmp_vals(arr1[i], arr2[i], eps);


    return ret;
}
/* ------------------------------------------------------------------------- */
