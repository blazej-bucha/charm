/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../prec.h"
#include "shc_reset_coeffs.h"
#include "shc_read_line.h"
#include "shc_read_gfc.h"
#include "../misc/misc_str2ul.h"
#include "../misc/misc_str2real.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






void CHARM(shc_read_gfc)(const char *pathname, unsigned long nmax,
                         CHARM(shc) *shcs, CHARM(err) *err)
{
    /* Open "pathname" to read */
    /* --------------------------------------------------------------------- */
    FILE *fptr = fopen(pathname, "r");
    if (fptr == NULL)
    {
        char msg[CHARM_ERR_MAX_MSG];
        sprintf(msg, "Couldn't open \"%s\".", pathname);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFILEIO, msg);
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Initial declarations */
    /* --------------------------------------------------------------------- */
    char line[SHC_READ_GFC_NLINE];
    char key_str[SHC_READ_GFC_NSTR];
    char val_str[SHC_READ_GFC_NSTR];
    unsigned long nmax_file;
    /* --------------------------------------------------------------------- */






    /* Read the header of the "gfc" file */
    /* --------------------------------------------------------------------- */
    _Bool nmax_found = 0;
    _Bool mu_found   = 0;
    _Bool r_found    = 0;
    _Bool eoh_found  = 0;
    int ret;


    do
    {

        /* Read a line of the file */
        if (fgets(line, SHC_READ_GFC_NLINE, fptr) == NULL)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                           CHARM_EFILEIO, "Couldn't find the \"end_of_head\" "
                                          "keyword or couldn't read from the "
                                          "input file.");
            goto EXIT;
        }


        /* Get the first two entries of "line" to test for the keywords of
         * "gfc" files */
        errno = 0;
        ret = sscanf(line, "%s %s", key_str, val_str);
        if (errno != 0)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                           CHARM_EFILEIO, "Couldn't read with \"sscanf\" from "
                                          "a file of the \"gfc\" file.");
            goto EXIT;
        }


        /* Is this the line specifying the end of header keyword? */
        if ((ret == 1 || ret == 2) && (strcmp(key_str, SHC_READ_GFC_EOH) == 0))
            eoh_found = 1;
        else
            eoh_found = 0;


        /* Is this the line specifying the maximum harmonic degree? */
        if ((ret == 2) && (strcmp(key_str, SHC_READ_GFC_NMAX) == 0))
        {
            nmax_file = CHARM(misc_str2ul)(val_str, "Failed to convert "
                                           "the maximum harmonic degree to "
                                           "the \"unsigned long int\" data "
                                           "format.", err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto EXIT;
            }
            nmax_found = 1;
        }


        /* Is this the line specifying the geocentric gravitational constant?
         * */
        if ((ret == 2) && (strcmp(key_str, SHC_READ_GFC_GM) == 0))
        {
            shcs->mu = CHARM(misc_str2real)(val_str, "Failed to convert "
                                            "the geocentric gravitational "
                                            "constant to the "
                                            "\"REAL\" data format.", err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto EXIT;
            }
            mu_found = 1;
        }


        /* Is this the line specifying the radius constant? */
        if ((ret == 2) && (strcmp(key_str, SHC_READ_GFC_R) == 0))
        {
            shcs->r = CHARM(misc_str2real)(val_str, "Failed to convert "
                                           "the radius constant to the "
                                           "\"REAL\" data format.", err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto EXIT;
            }
            r_found = 1;
        }


        /* Is this the line specifying the norm parameter? */
        if ((ret == 2) && (strcmp(key_str, SHC_READ_GFC_NORM) == 0))
        {
            if (strncmp(val_str, SHC_READ_GFC_NORM_FULL,
                        strlen(SHC_READ_GFC_NORM_FULL)) != 0)
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFILEIO, "Unsupported norm of spherical "
                               "harmonics.");
                goto EXIT;
            }


            /* "norm" is an optional keyword, so no "norm_found" here */
        }

    } while (!eoh_found);


    if (!nmax_found)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFILEIO, "Couldn't read the \"max_degree\" "
                       "keyword before reaching \"end_of_head\".");
        goto EXIT;
    }


    if (!mu_found)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFILEIO, "Couldn't read the "
                       "\"earth_gravity_constant\" keyword before reaching "
                       "\"end_of_head\".");
        goto EXIT;
    }


    if (!r_found)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFILEIO, "Couldn't read the "
                       "\"radius\" keyword before reaching "
                       "\"end_of_head\".");
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






    /* Read the table of spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    /* At first, reset all coefficients in "shcs" to zero. */
    CHARM(shc_reset_coeffs)(shcs);


    unsigned long n, m;
    REAL cnm, snm;
    while (CHARM(shc_read_line)(fptr, &n, &m, &cnm, &snm, SHC_READ_LINE_GFC,
                                err) != -1)
    {
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto EXIT;
        }


        if (n > nmax)
            continue;


        shcs->c[m][n - m] = cnm;
        shcs->s[m][n - m] = snm;
    }
    /* --------------------------------------------------------------------- */





EXIT:
    fclose(fptr);
    return;
}
