/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../prec.h"
#include "shc_read_line.h"
#include "../misc/misc_str2ul.h"
#include "../misc/misc_str2real.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






/* Reads and processes a single line of a table of spherical harmonic
 * coefficients from a text file or a gfc file.  Returns "-1" if EOF was
 * returned from "fgets" and "0" otherwise.  This function is not a part of
 * API. */
int CHARM(shc_read_line)(FILE *stream,
                         unsigned long *n,
                         unsigned long *m,
                         REAL *cnm,
                         REAL *snm,
                         int format,
                         CHARM(err) *err)
{
    /* Char arrays to store the strings to be loaded from the input file */
    /* --------------------------------------------------------------------- */
    char line[SHC_READ_LINE_NLINE];
    char n_str[SHC_READ_LINE_NSTR];
    char m_str[SHC_READ_LINE_NSTR];
    char cnm_str[SHC_READ_LINE_NSTR];
    char snm_str[SHC_READ_LINE_NSTR];
    /* --------------------------------------------------------------------- */






    /* Get the line of the input file */
    /* --------------------------------------------------------------------- */
    if (fgets(line, SHC_READ_LINE_NLINE, stream) == NULL)
        return -1;
    /* --------------------------------------------------------------------- */






    /* Extract strings from "line" and convert them to numerical data formats
     * */
    /* --------------------------------------------------------------------- */
    int num_entries;
    if (format == SHC_READ_LINE_TBL)
    {

        num_entries = sscanf(line, "%s %s %s %s",
                             n_str, m_str, cnm_str, snm_str);
        if ((num_entries != 3) && (num_entries != 4))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                           "Not enough entries in the coefficients "
                           "table line.");
            return 0;
        }

    }
    else if (format == SHC_READ_LINE_GFC)
    {

        char gfc_str[SHC_READ_LINE_NSTR];


        num_entries = sscanf(line, "%s %s %s %s %s",
                             gfc_str, n_str, m_str, cnm_str, snm_str);
        if (num_entries == -1)
            /* Probably an empty line which is valid in "gfc" files */
            return 0;


        if ((strcmp(gfc_str, SHC_READ_LINE_GFCT_KEYWORD) == 0) ||
            (strcmp(gfc_str, SHC_READ_LINE_TRND_KEYWORD) == 0) ||
            (strcmp(gfc_str, SHC_READ_LINE_ASIN_KEYWORD) == 0) ||
            (strcmp(gfc_str, SHC_READ_LINE_ACOS_KEYWORD) == 0))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                           "The \"gfc\" keyword is required in all "
                           "coefficients table lines (the data section).");
            return 0;
        }


        /* Any line starting with any text other keyword than
         * "SHC_READ_LINE_GFC_KEYWORD" is a comment, so should be skipped */
        if (strcmp(gfc_str, SHC_READ_LINE_GFC_KEYWORD) != 0)
            return 0;


        if ((num_entries != 4) && (num_entries != 5))
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                           "Not enough entries in the coefficients "
                           "table line (the data section).");
            return 0;
        }

    }
    else
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Wrong value of the \"format\" input parameter.");
        return 0;
    }


    *n = CHARM(misc_str2ul)(n_str, "Failed to convert harmonic degree "
                            "to the \"unsigned long int\" data format.", err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return 0;
    }


    *m = CHARM(misc_str2ul)(m_str, "Failed to convert harmonic order "
                            "to the \"unsigned long int\" data format.", err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return 0;
    }


    *cnm = CHARM(misc_str2real)(cnm_str, "Failed to convert the \"cnm\" "
                                "coefficient to the \"REAL\" "
                                "data format.", err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return 0;
    }


    /* Some tables of spherical harmonic coefficients omit the "sn0"
     * coefficients, as they do not exist */
    if (((format == SHC_READ_LINE_TBL) && (num_entries == 3)) ||
        ((format == SHC_READ_LINE_GFC) && (num_entries == 4)))
    {
        if (*m == 0)
            *snm = PREC(0.0);
        else
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                           "Wrong number of entries in the coefficients "
                           "table line.");


        return 0;
    }


    *snm = CHARM(misc_str2real)(snm_str, "Failed to convert the \"snm\" "
                                "coefficient to the \"REAL\" "
                                "data format.", err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return 0;
    }
    /* --------------------------------------------------------------------- */






    return 0;
}
