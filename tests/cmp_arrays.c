/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "cmp_vals.h"
#include "cmp_arrays.h"
/* ------------------------------------------------------------------------- */






/* Compares elements of two arrays up to some threshold */
/* ------------------------------------------------------------------------- */
long int cmp_arrays(REAL *arr1, REAL *arr2, size_t n, REAL eps)
{
    long int ret = 0;


    for (size_t i = 0; i < n; i++)
        ret += cmp_vals(arr1[i], arr2[i], eps);


    return ret;
}
/* ------------------------------------------------------------------------- */
