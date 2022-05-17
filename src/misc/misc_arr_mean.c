/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(misc_arr_mean)(const REAL *a, size_t na, CHARM(err) *err)
{
    if (na < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "The array size cannot be smaller than \"1\".");


        return (ADDP(0.0) / ADDP(0.0));
    }


    REAL mv = a[0]; /* Initialization */
    for (size_t i = 1; i < na; i++)
        mv += a[i];


    mv /= (REAL)na;


    return mv;
}
