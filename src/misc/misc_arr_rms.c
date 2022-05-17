/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(misc_arr_rms)(const REAL *a, size_t na, CHARM(err) *err)
{
    if (na < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The array size cannot be smaller than \"1\".");


        return (ADDP(0.0) / ADDP(0.0));
    }


    /* Initialization */
    REAL rms = a[0] * a[0];


    for (size_t i = 1; i < na; i++)
        rms += a[i] * a[i];


    rms /= (REAL)na;
    rms = SQRT(rms);


    return rms;
}
