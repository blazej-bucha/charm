/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "shc_reset_coeffs.h"
/* ------------------------------------------------------------------------- */






/* Function prototypes */
/* ------------------------------------------------------------------------- */
static int read_cnmsnm(FILE *, unsigned long, int, CHARM(shc) *);
/* ------------------------------------------------------------------------- */






void CHARM(shc_read_bin)(const char *pathname, unsigned long nmax,
                         CHARM(shc) *shcs, CHARM(err) *err)
{
    /* Open "pathname" to read */
    /* ===================================================================== */
    FILE *fptr = fopen(pathname, "rb");
    if (fptr == NULL)
    {
        char msg[CHARM_ERR_MAX_MSG];
        sprintf(msg, "Couldn't open \"%s\".", pathname);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFILEIO, msg);
        return;
    }
    /* ===================================================================== */






    /* Read the maximum harmonic degree, the scaling parameter and the radius
     * of the reference sphere */
    /* ===================================================================== */
    /* Variable to check whether "fread" was successful */
    int err_tmp;


    unsigned long nmax_file;
    err_tmp = fread(&nmax_file, sizeof(unsigned long), 1, fptr);
    if (err_tmp < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to read the maximum harmonic degree.");
        goto EXIT;
    }


    REAL mu;
    err_tmp = fread(&mu, sizeof(REAL), 1, fptr);
    if (err_tmp < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to read the scaling parameter.");
        goto EXIT;
    }


    REAL r;
    err_tmp = fread(&r, sizeof(REAL), 1, fptr);
    if (err_tmp < 1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to read the radius of the reference sphere.");
        goto EXIT;
    }
    /* ===================================================================== */






    /* Check maximum harmonic degrees */
    /* ===================================================================== */
    if (nmax > nmax_file)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFUNCARG,
                       "Not enough coefficients in the input file "
                       "for the maximum harmonic degree \"nmax\".");
        goto EXIT;
    }


    if (shcs->nmax < nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFUNCARG,
                       "Too low maximum degree \"shcs->nmax\" to read "
                       "coefficients up to degree \"nmax\".");
        goto EXIT;
    }
    /* ===================================================================== */






    /* ===================================================================== */
    /* Save the scaling constant and the radius of the reference sphere to the
     * "shcs" structure.  Note that we do not touch "shcs->nmax", "shcs->nc"
     * and "shcs->ns". */
    shcs->mu = mu;
    shcs->r  = r;


    /* Reset all coefficients in "shcs" to zero */
    CHARM(shc_reset_coeffs)(shcs);
    /* ===================================================================== */






    /* Read the "shcs->c" coefficients */
    /* ===================================================================== */
    err_tmp = read_cnmsnm(fptr, nmax_file, 0, shcs);
    if (err_tmp != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to read the \"C\" coefficients.");
        goto EXIT;
    }
    /* ===================================================================== */






    /* Read the "shcs->s" coefficients */
    /* ===================================================================== */
    err_tmp = read_cnmsnm(fptr, nmax_file, 1, shcs);
    if (err_tmp != 0)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to read the \"S\" coefficients.");
        goto EXIT;
    }
    /* ===================================================================== */






EXIT:
    fclose(fptr);
    return;
}






/* Just a small function to read "Cnm" and "Snm" coefficients from the binary
 * file.  Hopefully, no detailed documentation is needed. */
static int read_cnmsnm(FILE *fptr, unsigned long nmax_file, int cnmsnm,
                       CHARM(shc) *shcs)
{
    int err = 0;


    /* "fseek" needs "move_ptr" to be of "long" data type, so we do the
     * conversion here explicitly, rather than implicitly when calling
     * "fseek". */
    long move_ptr = (long)((nmax_file - shcs->nmax) * sizeof(REAL));
    _Bool move_ptr_bool = (shcs->nmax < nmax_file);


    /* Loop over the harmonic orders.  Note that we have to iterate up to
     * degree "nmax_file" in order to properly read coefficients if "shcs->nmax
     * < nmax_file". */
    for (unsigned long m = 0; m <= nmax_file; m++)
    {
        if (m <= shcs->nmax)
        {
            if (cnmsnm == 0)
                /* We are reading the "C" coefficients */
                err = fread(shcs->c[m] + m, sizeof(REAL), nmax_file + 1 - m,
                            fptr);
            else if (cnmsnm == 1)
                /* We are reading the "S" coefficients */
                err = fread(shcs->s[m] + m, sizeof(REAL), nmax_file + 1 - m,
                            fptr);
            else
                return 1;


            if (err < 1)
                return 2;


            if (move_ptr_bool)
            {
                /* We are reading only a subset of coefficients from file.
                 * More specifically, instead of reading up to "nmax_file", we
                 * are reading only up to degree "shcs->nmax".  Therefore, we
                 * have to move the "fptr" pointer properly, so that we can
                 * read the correct data in the next iteration. */
                err = fseek(fptr, move_ptr, SEEK_CUR);
                if (err != 0)
                    return 3;
            }
        }
        else
        {
            /* We are not interested in coefficients beyond degree
             * "shcs->nmax", so let's skip them. */
            err = fseek(fptr, (long)((nmax_file + 1 - m) * sizeof(REAL)),
                        SEEK_CUR);
            if (err != 0)
                return 3;
        }
    }


    return 0;
}
