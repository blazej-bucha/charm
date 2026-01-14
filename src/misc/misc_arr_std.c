/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "misc_arr_mean.h"
#include "misc_arr_std.h"
#include "misc_nan.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(misc_arr_std)(const REAL *a, size_t na, CHARM(err) *err)
{
    if (na < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The array size cannot be smaller than \"1\".");
        return CHARM_NAN;
    }


    /* Find the mean value of "a" */
    REAL mean = CHARM(misc_arr_mean)(a, na, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return CHARM_NAN;
    }


    /* Initialization */
    REAL res  = a[0] - mean;
    REAL std  = res * res;


    for (size_t i = 1; i < na; i++)
    {
        res = a[i] - mean;
        std += res * res;
    }


    std /= (REAL)na;
    std = SQRT(std);


    return std;
}
