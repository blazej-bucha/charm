/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../prec.h"
#include "shc_reset_coeffs.h"
#include "shc_read_mtdt.h"
#include "shc_read_nmax_only.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../misc/misc_str2ul.h"
#include "../misc/misc_str2real.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants */
/* ------------------------------------------------------------------------- */
/* Size of char arrays to store values of values loaded from the "dov" file */
#undef SHC_READ_DOV_NSTR
#define SHC_READ_DOV_NSTR (128)


/* Size of the char array to store a single line of the "dov" file */
#undef SHC_READ_DOV_NLINE
#define SHC_READ_DOV_NLINE (2048)
/* ------------------------------------------------------------------------- */






unsigned long CHARM(shc_read_dov)(const char *pathname,
                                  unsigned long nmax,
                                  CHARM(shc) *shcs,
                                  CHARM(err) *err)
{
    /* Open "pathname" to read */
    /* --------------------------------------------------------------------- */
    FILE *fptr = fopen(pathname, "r");
    if (fptr == NULL)
    {
        char msg[CHARM_ERR_MAX_MSG];
        sprintf(msg, "Couldn't open \"%s\".", pathname);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO, msg);
        return CHARM_SHC_NMAX_ERROR;
    }
    /* --------------------------------------------------------------------- */






    /* Read the metadata of spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    unsigned long nmax_file = CHARM_SHC_NMAX_ERROR;
    REAL r_file, mu_file;


    CHARM(shc_read_mtdt)(fptr, &nmax_file, &mu_file, &r_file, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        goto EXIT;
    }


    if (CHARM(shc_read_nmax_only)(nmax, shcs))
        goto EXIT;


    shcs->mu = mu_file;
    shcs->r  = r_file;
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
    char line[SHC_READ_DOV_NLINE];
    char n_str[SHC_READ_DOV_NSTR];
    char m_str[SHC_READ_DOV_NSTR];
    char coeff_str[SHC_READ_DOV_NSTR];


    /* At first, reset all coefficients in "shcs" to zero. */
    CHARM(shc_reset_coeffs)(shcs);


    unsigned long n, m;
    /* Pointer to the minus sign in the order entry of the "dov" file.
     * "m_sign" is "NULL" if the order is positive. */
    char *m_sign;
    int num_entries;
    REAL coeff;
    while (fgets(line, SHC_READ_DOV_NLINE, fptr) != NULL)
    {
        errno = 0;
        num_entries = sscanf(line, "%s %s %s",
                             n_str, m_str, coeff_str);
        if (errno)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                           CHARM_EFILEIO, "Couldn't read with \"sscanf\" from "
                                          "the \"dov\" file.");
            goto EXIT;
        }


        if (num_entries != 3)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                           "Not enough entries in the coefficients "
                           "table line.");
            goto EXIT;
        }


        n = CHARM(misc_str2ul)(n_str, "Failed to convert harmonic degree "
                               "to the \"unsigned long int\" data format.",
                               err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto EXIT;
        }


        if (n > nmax)
            continue;


        /* In "dov" files, the harmonic order we need to read may be negative.
         * To avoid creating a new function "misc_str2ll" ("long long"), we
         * check "m_str" for the presence of the first minus sign.  If any, we
         * change the sign to plus, so that "misc_str2ul" can be used to read
         * the modified "m_str" as "unsigned long".  The value of "m_sign"
         * pointer is later used to decide whether "coeff" should be stored in
         * "shcs->cnm" or "shcs->snm".
         *
         * It is sufficient to check "m_str" for the presence of the first
         * minus sign.  Otherwise, "misc_str2ul" will throw an error anyway. */
        m_sign = strchr(m_str, '-');
        if (m_sign != NULL)
            *m_sign = '+';
        m = CHARM(misc_str2ul)(m_str, "Failed to convert harmonic order "
                               "to the \"unsigned long int\" data format.",
                               err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto EXIT;
        }


        coeff = CHARM(misc_str2real)(coeff_str, "Failed to convert the "
                                     "spherical harmonic coefficient to the "
                                     "\"REAL\" data format.", err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto EXIT;
        }


        if (m_sign == NULL)
            shcs->c[m][n - m] = coeff;
        else
            shcs->s[m][n - m] = coeff;
    }
    /* --------------------------------------------------------------------- */






EXIT:
    fclose(fptr);
    return nmax_file;
}
