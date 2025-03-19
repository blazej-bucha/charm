/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "../prec.h"
#include "shc_reset_coeffs.h"
#include "shc_read_mtdt.h"
#include "shc_read_nmax_only.h"
#include "shc_check_distribution.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_check_distribution.h"
#include "../misc/misc_str2ul.h"
#include "../misc/misc_str2real.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants */
/* ------------------------------------------------------------------------- */
/* Size of char arrays to store values of values loaded from the "tbl" file */
#undef SHC_READ_TBL_NSTR
#define SHC_READ_TBL_NSTR (128)


/* Size of the char array to store a single line of the "tbl" file */
#undef SHC_READ_TBL_NLINE
#define SHC_READ_TBL_NLINE (2048)
/* ------------------------------------------------------------------------- */






unsigned long CHARM(shc_read_tbl)(const char *pathname,
                                  unsigned long nmax,
                                  CHARM(shc) *shcs,
                                  CHARM(err) *err)
{
    /* --------------------------------------------------------------------- */
    CHARM(err_check_distribution)(err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return CHARM_SHC_NMAX_ERROR;
    }


    if (!CHARM(shc_read_nmax_only)(nmax, shcs))
    {
        CHARM(shc_check_distribution)(shcs, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return CHARM_SHC_NMAX_ERROR;
        }
    }
    /* --------------------------------------------------------------------- */






    /* Open "pathname" to read */
    /* --------------------------------------------------------------------- */
    FILE *fptr = fopen(pathname, "r");
    if (fptr == NULL)
    {
        char msg[CHARM_ERR_MAX_MSG];
        snprintf(msg, CHARM_ERR_MAX_MSG, "Couldn't open \"%s\".", pathname);
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
    char line[SHC_READ_TBL_NLINE];
    char n_str[SHC_READ_TBL_NSTR];
    char m_str[SHC_READ_TBL_NSTR];
    char cnm_str[SHC_READ_TBL_NSTR];
    char snm_str[SHC_READ_TBL_NSTR];


    /* At first, reset all coefficients in "shcs" to zero. */
    CHARM(shc_reset_coeffs)(shcs);


    unsigned long n, m;
    int num_entries;
    REAL cnm, snm;
    while (fgets(line, SHC_READ_TBL_NLINE, fptr) != NULL)
    {
        errno = 0;
        num_entries = sscanf(line, "%s %s %s %s",
                             n_str, m_str, cnm_str, snm_str);
        if (errno)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                           CHARM_EFILEIO, "Couldn't read with \"sscanf\" from "
                                          "the \"tbl\" file.");
            goto EXIT;
        }


        if ((num_entries != 3) && (num_entries != 4))
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


        m = CHARM(misc_str2ul)(m_str, "Failed to convert harmonic order "
                               "to the \"unsigned long int\" data format.",
                               err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto EXIT;
        }


        cnm = CHARM(misc_str2real)(cnm_str, "Failed to convert the \"cnm\" "
                                   "coefficient to the \"REAL\" "
                                   "data format.", err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            goto EXIT;
        }


        /* Some tables of spherical harmonic coefficients omit the "sn0"
         * coefficients, as they do not exist */
        if (num_entries == 3)
        {
            if (m == 0)
                snm = PREC(0.0);
            else
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFILEIO,
                               "Wrong number of entries in the coefficients "
                               "table line.");
                goto EXIT;
            }
        }
        else
        {
            snm = CHARM(misc_str2real)(snm_str, "Failed to convert the "
                                       "\"snm\" coefficient to the \"REAL\" "
                                       "data format.", err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto EXIT;
            }
        }


        shcs->c[m][n - m] = cnm;
        shcs->s[m][n - m] = snm;
    }
    /* --------------------------------------------------------------------- */






EXIT:
    fclose(fptr);
    return nmax_file;
}
