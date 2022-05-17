/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../misc/misc_fprintf_real.h"
#include "shc_write_mtdt.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_write_tbl)(const CHARM(shc) *shcs, unsigned long nmax,
                          const char *format, int order, FILE *stream,
                          CHARM(err) *err)
{
    /* Check maximum harmonic degree */
    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Not enough coefficients in \"shcs\" to write "
                       "up to degree \"nmax\".");
        return;
    }


    /* Write the metadata */
    CHARM(shc_write_mtdt)(nmax, shcs->mu, shcs->r, format, stream, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    /* Write the spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    if (order == CHARM_SHC_WRITE_TBL_N)
    {
        for (unsigned long m = 0; m <= nmax; m++)
        {
            for (unsigned long n = m; n <= nmax; n++)
            {
                if ((fprintf(stream, "%lu ", n) < 1) ||
                    (fprintf(stream, "%lu ", m) < 1) ||
                    (CHARM(misc_fprintf_real)(stream, format,
                                              shcs->c[m][n - m]) < 1) ||
                    (fprintf(stream, " ") < 1) ||
                    (CHARM(misc_fprintf_real)(stream, format,
                                              shcs->s[m][n - m]) < 1) ||
                    (fprintf(stream, "\n") < 1))
                {
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EFILEIO,
                                   "Failed writing to the output file.");
                    return;
                }
            }
        }
    }
    else if (order == CHARM_SHC_WRITE_TBL_M)
    {
        for (unsigned long n = 0; n <= nmax; n++)
        {
            for (unsigned long m = 0; m <= n; m++)
            {
                if ((fprintf(stream, "%lu ", n) < 1) ||
                    (fprintf(stream, "%lu ", m) < 1) ||
                    (CHARM(misc_fprintf_real)(stream, format,
                                              shcs->c[m][n - m]) < 1) ||
                    (fprintf(stream, " ") < 1) ||
                    (CHARM(misc_fprintf_real)(stream, format,
                                              shcs->s[m][n - m]) < 1) ||
                    (fprintf(stream, "\n") < 1))
                {
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EFILEIO,
                                   "Failed writing to the output file.");
                    return;
                }
            }
        }
    }
    else
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFUNCARG, "Failed writing to the output file.");
        return;
    }
    /* --------------------------------------------------------------------- */


    return;
}
