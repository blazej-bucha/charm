/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shc_check_distribution.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_check_distribution.h"
/* ------------------------------------------------------------------------- */






/* Function prototypes */
/* ------------------------------------------------------------------------- */
static int write_cnmsnm(const CHARM(shc) *,
                        unsigned long,
                        _Bool,
                        FILE *);
/* ------------------------------------------------------------------------- */






void CHARM(shc_write_bin)(const CHARM(shc) *shcs,
                          unsigned long nmax,
                          const char *pathname,
                          CHARM(err) *err)
{
    /* ===================================================================== */
    CHARM(err_check_distribution)(err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    CHARM(shc_check_distribution)(shcs, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }
    /* ===================================================================== */






    /* Open "pathname" to write */
    /* ===================================================================== */
    FILE *fptr = fopen(pathname, "wb");
    if (fptr == NULL)
    {
        char msg[CHARM_ERR_MAX_MSG];
        sprintf(msg, "Couldn't create \"%s\".", pathname);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFILEIO, msg);
        return;
    }
    /* ===================================================================== */






    /* Check maximum harmonic degree */
    /* ===================================================================== */
    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Not enough coefficients in \"shcs\" to write "
                       "up to degree \"nmax\".");
        goto EXIT;
    }
    /* ===================================================================== */






    /* Write the maximum harmonic degree, the scaling parameter and the radius
     * of the reference sphere */
    /* ===================================================================== */
    /* The maximum harmonic degree.  Note that written is the user-defined
     * maximum harmonic degree "nmax" (not "shcs->nmax"), as this is what we
     * are asked to do by the user. */
    if (fwrite(&nmax, sizeof(unsigned long), 1, fptr) < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the maximum harmonic degree.");
        goto EXIT;
    }


    /* The scaling parameter */
    if (fwrite(&(shcs->mu), sizeof(REAL), 1, fptr) < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the scaling parameter.");
        goto EXIT;
    }


    /* The scaling parameter */
    if (fwrite(&(shcs->r), sizeof(REAL), 1, fptr) < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the radius of the reference sphere.");
        goto EXIT;
    }
    /* ===================================================================== */






    /* Write the "shcs->c" coefficients */
    /* ===================================================================== */
    if (write_cnmsnm(shcs, nmax, 0, fptr) != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the \"C\" coefficients).");
        goto EXIT;
    }
    /* ===================================================================== */






    /* Write the "shcs->s" coefficients */
    /* ===================================================================== */
    if (write_cnmsnm(shcs, nmax, 1, fptr) != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the \"S\" coefficients).");
        goto EXIT;
    }
    /* ===================================================================== */






EXIT:
    fclose(fptr);
    return;
}






/* Just a small function to write "Cnm" and "Snm" coefficients to the binary
 * file.  Hopefully, no detailed documentation is needed. */
static int write_cnmsnm(const CHARM(shc) *shcs,
                        unsigned long nmax,
                        _Bool cnmsnm,
                        FILE *fptr)
{
    /* Loop over the harmonic orders */
    for (unsigned long m = 0; m <= nmax; m++)
    {
        REAL *cs = (cnmsnm) ? shcs->s[m] : shcs->c[m];
        if (fwrite(cs, sizeof(REAL), nmax + 1 - m, fptr) < 1)
            return 1;
    }


    return 0;
}
