/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "shc_reset_coeffs.h"
#include "shc_read_mtdt.h"
#include "../misc/misc_str2real.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants */
/* ------------------------------------------------------------------------- */
#undef NCHARS
#define NCHARS (256)
/* ------------------------------------------------------------------------- */






void CHARM(shc_read_mtx)(const char *pathname, unsigned long nmax,
                         CHARM(shc) *shcs, CHARM(err) *err)
{
    /* Open "pathname" to read */
    /* ===================================================================== */
    FILE *fptr = fopen(pathname, "r");
    if (fptr == NULL)
    {
        char msg[CHARM_ERR_MAX_MSG];
        sprintf(msg, "Couldn't open \"%s\".", pathname);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFILEIO, msg);
        return;
    }
    /* ===================================================================== */






    /* Read file */
    /* ===================================================================== */
    /* A string to be loaded from the input file */
    char str[NCHARS];


    /* String to store the new line character */
    char nl[] = "\0";


    /* The number of input items successfully matched and assigned from the
     * "fscanf" function */
    int num_entries;


    /* An entry of REAL data type from "fptr" */
    REAL entry_d;






    /* Read the metadata of spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    unsigned long nmax_file;
    CHARM(shc_read_mtdt)(fptr, &nmax_file, &(shcs->mu), &(shcs->r), err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */







    /* Check maximum harmonic degrees */
    /* --------------------------------------------------------------------- */
    if (shcs->nmax < nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Too low maximum degree \"shcs->nmax\" to read "
                       "coefficients up to degree \"nmax\".");
        goto EXIT;
    }


    if (nmax_file < nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Too low maximum degree inside the input file to read "
                       "coefficients up to degree \"nmax\".");
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* Read the spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    /* At first, reset all coefficients in "shcs" to zero. */
    CHARM(shc_reset_coeffs)(shcs);


    for (unsigned long row = 0; row <= nmax; row++)
    {
        /* Loop over the columns of the matrix */
        for (unsigned long col = 0; col <= nmax; col++)
        {
            /* Read an entry from the text file and store it as a string */
            /* ------------------------------------------------------------- */
            num_entries = fscanf(fptr, "%s", str);


            /* Check for EOF */
            if (num_entries == EOF)
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFILEIO,
                               "Too few rows in the input file to read "
                               "spherical harmonic coefficients up to degree "
                               "\"nmax\".");
                goto EXIT;
            }


            if (num_entries < 1)
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFILEIO,
                               "Failed to read an entry from the input file.");
                goto EXIT;
            }
            /* ------------------------------------------------------------- */


            /* Read the new line character if any */
            /* ------------------------------------------------------------- */
            num_entries = fscanf(fptr, "%1[\n]", nl);
            if ((num_entries == 1) && (nl[0] == '\n') && (col < nmax))
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFILEIO,
                               "Too few columns to read spherical harmonic "
                               "coefficients up to degree \"nmax\".");
                goto EXIT;
            }
            /* ------------------------------------------------------------- */


            /* Convert the number in "str" from string to REAL and check the
             * conversion */
            /* ------------------------------------------------------------- */
            entry_d = CHARM(misc_str2real)(str, "Failed to convert an entry "
                                           "from the input file to the "
                                           "\"REAL\" data format.", err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto EXIT;
            }
            /* ------------------------------------------------------------- */


            /* Save the value of "entry" to the right place of either "shcs->c"
             * or "shcs->s" */
            /* ------------------------------------------------------------- */
            if (row >= col)
                shcs->c[col][row] = entry_d;
            else
                shcs->s[row + 1][col] = entry_d;
            /* ------------------------------------------------------------- */
        }


        /* If the input data file contains more than "nmax + 1" columns,
         * continue with reading the next line, since all sought "nmax + 1"
         * values from the "row"th line were read */
        if (nl[0] != '\n')
            fscanf(fptr, "%*[^\n]\n");


        /* Reset "nl" (it may contain the new line character found at the end
         * of a line) */
        nl[0] = '\0';
    }
    /* --------------------------------------------------------------------- */
    /* ===================================================================== */






EXIT:
    fclose(fptr);
    return;
}

