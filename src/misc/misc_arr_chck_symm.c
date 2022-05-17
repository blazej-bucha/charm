/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_arr_chck_symm)(const REAL *arr, size_t size, REAL center,
                              REAL eps, CHARM(err) *err)
{
    if (size < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The array size cannot be smaller than \"1\".");


        return -9999;
    }


    /* Get the parity of "size"; "0" if even, "1" if odd. */
    int parity = size % 2;


    /* Check the first condition, that is, symmetry of "arr" */
    size_t middle = size / 2;
    for (size_t i = 0; i < middle; i++)
        if (!CHARM(misc_is_nearly_equal)(arr[i] - center,
                                         center - arr[size - 1 - i],
                                         eps))
            return 1;


    /* If "size" is odd, check whether the middle element of "arr" equals to
     * "center" */
    if (parity == 1)
        if (CHARM(misc_is_nearly_equal)(center, arr[middle], eps) == 0)
            return 2;


    /* "arr" passed all tests successfully */
    return 0;
}
