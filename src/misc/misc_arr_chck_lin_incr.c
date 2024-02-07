/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "misc_is_nearly_equal.h"
#include "misc_arr_min.h"
#include "misc_arr_chck_lin_incr.h"
/* ------------------------------------------------------------------------- */






int CHARM(misc_arr_chck_lin_incr)(const REAL *arr, size_t size,
                                  size_t first, size_t every_nth,
                                  REAL eps, CHARM(err) *err)
{
    /* Error checks */
    /* --------------------------------------------------------------------- */
    if (size < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The array size cannot be smaller than \"1\".");


        return -9999;
    }


    if (every_nth < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The \"every_nth\" value cannot be smaller than "
                       "\"1\".");


        return -9999;
    }


    if (first >= size)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"first\" must be smaller than \"size\".");


        return -9999;

    }


    if (every_nth >= size)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"every_nth\" must be smaller than \"size\".");


        return -9999;

    }
    /* --------------------------------------------------------------------- */


    /* Check the first condition */
    /* --------------------------------------------------------------------- */
    if (CHARM(misc_is_nearly_equal)(arr[first],
                                    CHARM(misc_arr_min)(arr + first,
                                                        size - first,
                                                        err),
                                    eps) == 0)
        return 1;
    /* --------------------------------------------------------------------- */


    /* Get the step of the array considering the "every_nth" value */
    /* --------------------------------------------------------------------- */
    REAL step;
    if (size > 1)
    {
        step = arr[first + every_nth] - arr[first];

        /* Check whether "step" is positive */
        if (step <= PREC(0.0))
            return 2;
    }
    else
        step = PREC(0.0);
    /* --------------------------------------------------------------------- */


    /* Check whether the array step is constant for the entire array. */
    /* --------------------------------------------------------------------- */
    for (size_t i = first + every_nth; i < size; i += every_nth)
        if (CHARM(misc_is_nearly_equal)(step, arr[i] - arr[i - every_nth],
                                        eps) == 0)
            return 2;
    /* --------------------------------------------------------------------- */


    return 0;
}
