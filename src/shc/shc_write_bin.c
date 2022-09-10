/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






/* Function prototypes */
/* ------------------------------------------------------------------------- */
static int write_cnmsnm(const CHARM(shc) *, unsigned long, int, FILE *);
/* ------------------------------------------------------------------------- */






void CHARM(shc_write_bin)(const CHARM(shc) *shcs, unsigned long nmax,
                          const char *pathname, CHARM(err) *err)
{
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
    /* Variable to check whether "fwrite" was successful */
    int err_tmp;


    /* The maximum harmonic degree.  Note that written is the user-defined
     * maximum harmonic degree "nmax" (not "shcs->nmax"), as this is what we
     * are asked to do by the user. */
    err_tmp = fwrite(&nmax, sizeof(unsigned long), 1, fptr);
    if (err_tmp < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the maximum harmonic degree.");
        goto EXIT;
    }


    /* The scaling parameter */
    err_tmp = fwrite(&(shcs->mu), sizeof(REAL), 1, fptr);
    if (err_tmp < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the scaling parameter.");
        goto EXIT;
    }


    /* The scaling parameter */
    err_tmp = fwrite(&(shcs->r), sizeof(REAL), 1, fptr);
    if (err_tmp < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the radius of the reference sphere.");
        goto EXIT;
    }
    /* ===================================================================== */






    /* Write the "shcs->c" coefficients */
    /* ===================================================================== */
    err_tmp = write_cnmsnm(shcs, nmax, 0, fptr);
    if (err_tmp != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the \"C\" coefficients).");
        goto EXIT;
    }
    /* ===================================================================== */






    /* Write the "shcs->s" coefficients */
    /* ===================================================================== */
    err_tmp = write_cnmsnm(shcs, nmax, 1, fptr);
    if (err_tmp != 0)
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
static int write_cnmsnm(const CHARM(shc) *shcs, unsigned long nmax, int cnmsnm,
                        FILE *fptr)
{
    int err = 0;


    /* Loop over the harmonic orders */
    for (unsigned long m = 0; m <= nmax; m++)
    {
        if (cnmsnm == 0)
            err = fwrite(shcs->c[m], sizeof(REAL), nmax + 1 - m, fptr);
        else if (cnmsnm == 1)
            err = fwrite(shcs->s[m], sizeof(REAL), nmax + 1 - m, fptr);
        else
            return 1;


        if (err < 1)
            return 2;
    }


    return 0;
}
