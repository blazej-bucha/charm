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






void CHARM(shc_write_tbl)(const CHARM(shc) *shcs,
                          unsigned long nmax,
                          const char *format,
                          int ordering,
                          const char *pathname,
                          CHARM(err) *err)
{
    /* Open "pathname" to write */
    /* --------------------------------------------------------------------- */
    FILE *fptr = fopen(pathname, "w");
    if (fptr == NULL)
    {
        char msg[CHARM_ERR_MAX_MSG];
        sprintf(msg, "Couldn't create \"%s\".", pathname);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFILEIO, msg);
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Check maximum harmonic degree */
    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Not enough coefficients in \"shcs\" to write "
                       "up to degree \"nmax\".");
        goto EXIT;
    }


    /* Write the metadata */
    CHARM(shc_write_mtdt)(nmax, shcs->mu, shcs->r, format, fptr, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }


    /* Write the spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    if (ordering == CHARM_SHC_WRITE_N)
    {
        for (unsigned long m = 0; m <= nmax; m++)
        {
            for (unsigned long n = m; n <= nmax; n++)
            {
                if ((fprintf(fptr, "%lu ", n) < 1) ||
                    (fprintf(fptr, "%lu ", m) < 1) ||
                    (CHARM(misc_fprintf_real)(fptr, format,
                                              shcs->c[m][n - m]) < 1) ||
                    (fprintf(fptr, " ") < 1) ||
                    (CHARM(misc_fprintf_real)(fptr, format,
                                              shcs->s[m][n - m]) < 1) ||
                    (fprintf(fptr, "\n") < 1))
                {
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EFILEIO,
                                   "Failed writing to the output file.");
                    goto EXIT;
                }
            }
        }
    }
    else if (ordering == CHARM_SHC_WRITE_M)
    {
        for (unsigned long n = 0; n <= nmax; n++)
        {
            for (unsigned long m = 0; m <= n; m++)
            {
                if ((fprintf(fptr, "%lu ", n) < 1) ||
                    (fprintf(fptr, "%lu ", m) < 1) ||
                    (CHARM(misc_fprintf_real)(fptr, format,
                                              shcs->c[m][n - m]) < 1) ||
                    (fprintf(fptr, " ") < 1) ||
                    (CHARM(misc_fprintf_real)(fptr, format,
                                              shcs->s[m][n - m]) < 1) ||
                    (fprintf(fptr, "\n") < 1))
                {
                    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                                   CHARM_EFILEIO,
                                   "Failed writing to the output file.");
                    goto EXIT;
                }
            }
        }
    }
    else
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFUNCARG, "Unsupported value of \"ordering\".");
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */


EXIT:
    fclose(fptr);
    return;
}
