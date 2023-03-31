/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/misc/misc_is_nearly_equal.h"
#include "parameters.h"
#include "cmp_vals.h"
/* ------------------------------------------------------------------------- */






/* Compares elements of two arrays up to some threshold */
/* ------------------------------------------------------------------------- */
int cmp_arrays(REAL *arr1, REAL *arr2, size_t n, REAL eps)
{
    int ret = 0;


    for (size_t i = 0; i < n; i++)
        ret += cmp_vals(arr1[i], arr2[i], eps);


    return ret;
}
/* ------------------------------------------------------------------------- */
