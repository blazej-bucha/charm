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






void CHARM(shc_write_mtx)(const CHARM(shc) *shcs, unsigned long nmax,
                          const char *format, const char *pathname,
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
    /* A particular coefficient to be written to the text file, either
     * "shcs->c" or "shcs->s" */
    REAL coeff;


    /* Loop over the rows of the output matrix */
    for (unsigned long row = 0; row <= nmax; row++)
    {
        /* Loop over the columns of the output matrix */
        for (unsigned long col = 0; col <= nmax; col++)
        {
            if (row >= col)
                coeff = shcs->c[col][row];
            else
                coeff = shcs->s[row + 1][col];


            if (CHARM(misc_fprintf_real)(fptr, format, coeff) < 1)
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFILEIO,
                               "Failed to write to the output file.");
                goto EXIT;
            }


            if (col < nmax)
                fprintf(fptr, " ");
        }


        if (fprintf(fptr, "\n") < 1)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                           CHARM_EFILEIO,
                           "Failed to write to the output file.");
            goto EXIT;
        }
    }
    /* --------------------------------------------------------------------- */


EXIT:
    fclose(fptr);
    return;
}
