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
                          FILE *stream, CHARM(err) *err)
{
    /* Check maximum harmonic degree */
    /* ===================================================================== */
    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Not enough coefficients in \"shcs\" to write "
                       "up to degree \"nmax\".");
        return;
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
    err_tmp = fwrite(&nmax, sizeof(unsigned long), 1, stream);
    if (err_tmp < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the maximum harmonic degree.");
        return;
    }


    /* The scaling parameter */
    err_tmp = fwrite(&(shcs->mu), sizeof(REAL), 1, stream);
    if (err_tmp < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the scaling parameter.");
        return;
    }


    /* The scaling parameter */
    err_tmp = fwrite(&(shcs->r), sizeof(REAL), 1, stream);
    if (err_tmp < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the radius of the reference sphere.");
        return;
    }
    /* ===================================================================== */






    /* Write the "shcs->c" coefficients */
    /* ===================================================================== */
    err_tmp = write_cnmsnm(shcs, nmax, 0, stream);
    if (err_tmp != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the \"C\" coefficients).");
        return;
    }
    /* ===================================================================== */






    /* Write the "shcs->s" coefficients */
    /* ===================================================================== */
    err_tmp = write_cnmsnm(shcs, nmax, 1, stream);
    if (err_tmp != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to write the \"S\" coefficients).");
        return;
    }
    /* ===================================================================== */






    return;
}






/* Just a small function to write "Cnm" and "Snm" coefficients to the binary
 * file.  Hopefully, no detailed documentation is needed. */
static int write_cnmsnm(const CHARM(shc) *shcs, unsigned long nmax, int cnmsnm,
                        FILE *stream)
{
    int err = 0;


    /* Loop over the harmonic orders */
    for (unsigned long m = 0; m <= nmax; m++)
    {
        if (cnmsnm == 0)
            err = fwrite(shcs->c[m], sizeof(REAL), nmax + 1 - m, stream);
        else if (cnmsnm == 1)
            err = fwrite(shcs->s[m], sizeof(REAL), nmax + 1 - m, stream);
        else
            return 1;


        if (err < 1)
            return 2;
    }


    return 0;
}
