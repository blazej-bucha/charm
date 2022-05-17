/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../misc/misc_fprintf_real.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_write_mtdt)(unsigned long nmax, REAL mu, REAL r,
                           const char *format, FILE *stream, CHARM(err) *err)
{
    /* Write the maximum harmonic degree */
    if (fprintf(stream, "%lu ", nmax) < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the maximum harmonic degree.");
        return;
    }


    /* Write the scaling parameter */
    if (CHARM(misc_fprintf_real)(stream, format, mu) < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the scaling parameter.");
        return;
    }
    fprintf(stream, " ");


    /* Write the radius of the reference sphere */
    if (CHARM(misc_fprintf_real)(stream, format, r) < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the radius of the reference sphere.");
        return;
    }
    fprintf(stream, "\n");


    return;
}
