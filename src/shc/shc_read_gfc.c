/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include "../prec.h"
#include "shc_reset_coeffs.h"
#include "shc_read_gfc.h"
#include "shc_read_nmax_only.h"
#include "../misc/misc_str2ul.h"
#include "../misc/misc_str2real.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants */
/* ------------------------------------------------------------------------- */
#undef DAYS_PER_NONLEAP_YEAR
#define DAYS_PER_NONLEAP_YEAR (PREC(365.0))


#undef HOURS_PER_DAY
#define HOURS_PER_DAY (PREC(24.0))


#undef MINUTES_PER_HOUR
#define MINUTES_PER_HOUR (PREC(60.0))


#undef HOUR_MAX
#define HOUR_MAX (24)


#undef MINUTE_MAX
#define MINUTE_MAX (60)


/* Return this value if "epoch_fraction" fails */
#undef EPOCH_FRACTION_ERR
#define EPOCH_FRACTION_ERR (PREC(-9999.0))


/* Return this value if "epoch" in "epoch_fraction" is NULL */
#undef EPOCH_FRACTION_NULL
#define EPOCH_FRACTION_NULL (PREC(-8888.0))
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
/* Read harmonic degree/order from a line */
#define READ_DEG_ORD(val, str, degord, pathname)                              \
        val = CHARM(misc_str2ul)(str, "", err);                               \
        if (!CHARM(err_isempty)(err))                                         \
        {                                                                     \
            CHARM(err_reset)(err);                                            \
            sprintf(err_msg, "Failed to convert harmonic %s \"%s\" in \"%s\" "\
                             "to the \"unsigned long int\" data format.",     \
                    degord, str, pathname);                                   \
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,  \
                           err_msg);                                          \
            goto EXIT;                                                        \
        }


/* Read coefficients from a line */
#define READ_CNM_SNM(val, str, type, n, m, pathname)                          \
        val = CHARM(misc_str2real)(str, "", err);                             \
        if (!CHARM(err_isempty)(err))                                         \
        {                                                                     \
            CHARM(err_reset)(err);                                            \
            sprintf(err_msg, "Failed to convert the \"%s\" coefficient "      \
                             "\"%s\" of degree \"%lu\" and order \"%lu\" in " \
                             "\"%s\" to the \"REAL\" data format.",           \
                    type, str, n, m, pathname);                               \
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,  \
                           err_msg);                                          \
            goto EXIT;                                                        \
        }


/* Check the number of entries in a line */
#define WRONG_NUMBER_OF_ENTRIES                                               \
        remove_new_line_char(line);                                           \
        sprintf(err_msg, "The following line of \"%s\" has an "               \
                         "incorrect number of entries: \"%s\".",              \
                pathname, line);                                              \
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,                     \
                       CHARM_EFILEIO, err_msg);                               \
        goto EXIT;


/* In time variable models, checks whether degree and order of coefficients of
 * the types "SHC_READ_GFC_TRND", "SHC_READ_GFC_DOT", "SHC_READ_GFC_ASIN" and
 * "SHC_READ_GFC_ACOS" match their expected value. */
#define EXPECTED_DEGREE_ORDER(type, n, m, pathname, n_tmp, m_tmp, cnm, snm)   \
        if ((n_tmp != n) || (m_tmp != m))                                     \
        {                                                                     \
            sprintf(err_msg, "Expected in file \"%s\" coefficients \"%s\" of "\
                             "degree \"%lu\" and order \"%lu\", but "         \
                             "found instead coefficients \"%s\" and \"%s\" "  \
                             "of degree \"%lu\" and order \"%lu\".",          \
                    pathname, type, n, m, cnm, snm, n_tmp, m_tmp);            \
            CHARM(err_set)(err, __FILE__, __LINE__, __func__,                 \
                           CHARM_EFILEIO, err_msg);                           \
            goto EXIT;                                                        \
        }
/* ------------------------------------------------------------------------- */






/* Data types */
/* ------------------------------------------------------------------------- */
/* The data type to store the epoch as a fraction of year depends on the
 * precision, in which CHarm is being compiled.  The reason is that if
 * compiling in single precision, the range of the "float" data type may not be
 * sufficient to store the fraction of a year value.  In that case, double
 * precision is used.  When compiling CHarm in double or quadruple precision,
 * used is double or quadruple precision, respectively. */
#undef REAL_EPOCH
#if CHARM_FLOAT
#   define REAL_EPOCH double
#else
#   define REAL_EPOCH REAL
#endif
/* ------------------------------------------------------------------------- */






/* Functions */
/* ------------------------------------------------------------------------- */
static REAL_EPOCH epoch_fraction(const char *date,
                                 CHARM(err) *err)
{
    if (date == NULL)
        return EPOCH_FRACTION_NULL;


    char err_msg[CHARM_ERR_MAX_MSG];


    /* Let "date" to point to the first non-space character */
    while (isspace(date[0]))
        date++;


    /* BC dates are not allowed */
    if (date[0] == '-')
    {
        sprintf(err_msg, "Negative epoch \"%s\".", date);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto FAILURE;
    }


    size_t ndate = strlen(date);
    if (ndate < 8)
    {
        sprintf(err_msg, "Wrong format of the epoch string \"%s\".  "
                         "Supported formats are \"yyyyMMdd\" and "
                         "\"yyyyMMdd.hhmm\".", date);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto FAILURE;
    }


    /* The first value of the array is zero, because we have no zeroth month as
     * would be nice to have, given the indexing in C.  Then, the total number
     * of days per month (JAN, FEB, ..., DEC) follow. */
    static const unsigned long monthdays[13] = {0, 31, 28, 31, 30, 31, 30,
                                                31, 31, 30, 31, 30, 31};


    /* Zero and then leap days for a given month. */
    static const unsigned long leapdays[13] = {0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0};


    /* Zero and then the total number of days that have already passed before
     * the first of JAN, FEB, ..., DEC follow.  We do not consider here leap
     * years to allow the array to be static.  The leap years are treated
     * later. */
    static const REAL_EPOCH cumsumdays[13] = {PREC(0.0), PREC(0.0), PREC(31.0),
                                              PREC(59.0), PREC(90.0),
                                              PREC(120.0), PREC(151.0),
                                              PREC(181.0), PREC(212.0),
                                              PREC(243.0), PREC(273.0),
                                              PREC(304.0), PREC(334.0)};


    char yyyy[5], MM[3], dd[3], hh[3], mm[3];
    yyyy[4] = MM[2] = dd[2] = hh[2] = mm[2] = '\0';


    /* Split "date" to strings "yyyy", "MM", ... */
    memcpy(yyyy, date, 4);
    memcpy(MM,   date + 4, 2);
    memcpy(dd,   date + 6, 2);


    /* Convert "yyyy", etc to "REAL_EPOCH" or "unsigned long".  Instead of
     * "unsigned long", "int" or "unsigned int" could be used.  But given that
     * we already have a nice routine to convert a string to "unsigned long",
     * we use "unsigned int". */
    unsigned long year = CHARM(misc_str2ul)(yyyy, "", err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_reset)(err);
        sprintf(err_msg, "Failed to convert the year value \"yyyy = %s\" of "
                         "the epoch string \"%s\" to a floating point.", yyyy,
                         date);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto FAILURE;
    }


    /* Is this a leap year? */
    _Bool isleapyear = (!(year % 400) || (!(year % 4) && (year % 100))) ? 1
                                                                        : 0;


    unsigned long month = CHARM(misc_str2ul)(MM, "", err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_reset)(err);
        sprintf(err_msg, "Failed to convert the month value \"MM = %s\" of "
                         "the epoch string \"%s\" to a floating point.", MM,
                         date);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto FAILURE;
    }
    if ((month < 1) || (month > 12))
    {
        sprintf(err_msg, "Invalid value of month \"MM = %02lu\" in the "
                         "epoch string \"%s\".", month, date);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto FAILURE;
    }


    unsigned long day = CHARM(misc_str2ul)(dd, "", err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_reset)(err);
        sprintf(err_msg, "Failed to convert the day value \"dd = %s\" in the "
                         "epoch string \"%s\" to a floating point.", dd, date);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto FAILURE;
    }
    if ((day < 1) ||
        (day > (monthdays[month] + isleapyear * leapdays[month])))
    {
        sprintf(err_msg, "Invalid value of day \"day = %02lu\" in the "
                         "epoch string \"%s\".", day, date);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto FAILURE;
    }


    /* Check whether "date" is of "yyyyMMdd.hhmm" or "yyyyMMdd" structure */
    unsigned long hour, min;
    hour = min = 0;
    /* If the character after "yyyyMMdd" is not a terminating byte, it must
     * be ".", after which "hhmm" must follow. */
    if (ndate > 8)
    {
        if ((ndate != 13) || (date[8] != '.'))
        {
            sprintf(err_msg, "Wrong format of the epoch string \"%s\".  "
                             "Supported formats are \"yyyyMMdd\" and "
                             "\"yyyyMMdd.hhmm\".", date);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            goto FAILURE;
        }


        memcpy(hh, date + 9, 2);
        hour = CHARM(misc_str2ul)(hh, "", err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_reset)(err);
            sprintf(err_msg, "Failed to convert the value of hours "
                             "\"hh = %s\" in the epoch string string \"%s\" "
                             "to a floating point.", hh, date);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            goto FAILURE;
        }
        if (hour > HOUR_MAX)
        {
            sprintf(err_msg, "Invalid value of hour \"hh = %02lu\" in the "
                             "epoch string \"%s\".", hour, date);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            goto FAILURE;
        }


        memcpy(mm, date + 11, 2);
        min = CHARM(misc_str2ul)(mm, "", err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_reset)(err);
            sprintf(err_msg, "Failed to convert the value of minutes "
                             "\"mm = %s\" in the epoch string string \"%s\" "
                             "to a floating point.", mm, date);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            goto FAILURE;
        }
        if (min > MINUTE_MAX)
        {
            sprintf(err_msg, "Invalid value of minutes \"mm = %02lu\" in the "
                             "epoch string \"%s\".", min, date);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           err_msg);
            goto FAILURE;
        }
    }


    if ((hour == 24) && (min == 60))
    {
        sprintf(err_msg, "Invalid combination of hours and minutes "
                         "in the epoch string \"%s\".  If by \"hh = %02lu\" "
                         "and \"mm = %02lu\" you mean the beginning of the "
                         "next day, use either \"%04lu%02lu%02lu.0000\" or "
                         "\"%04lu%02lu%02lu.2400\".",
                         date, hour, min, year, month, day + 1, year, month,
                         day);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto FAILURE;
    }


    /* Day of year decreased by one */
    REAL_EPOCH doy = (REAL_EPOCH)day + cumsumdays[month] - PREC(1.0);
    if (isleapyear && (month > 2))
        doy += PREC(1.0);


    /* Total number of days in year */
    REAL_EPOCH ndays = DAYS_PER_NONLEAP_YEAR + isleapyear;


    return (REAL_EPOCH)year + (doy + ((REAL_EPOCH)hour +
                                      (REAL_EPOCH)min / MINUTES_PER_HOUR) /
                               HOURS_PER_DAY) / ndays;


FAILURE:
    return EPOCH_FRACTION_ERR;
}






/* Returns true if "t0 <= t < t1" and false otherwise. */
static _Bool is_within_time_period(REAL_EPOCH t,
                                   REAL_EPOCH t0,
                                   REAL_EPOCH t1)
{
    return (t >= t0) && (t < t1);
}






/* Computes the difference between two fractions of a year.  If needed (i.e.,
 * when compiling CHarm in single precision), the output difference is type
 * casted from "REAL_EPOCH" to "REAL". */
static REAL epoch_diff(REAL_EPOCH t, REAL_EPOCH t0)
{
    return t - t0;
}






/* Removes the first new line character in input string */
static void remove_new_line_char(char *str)
{
    char *replace = strchr(str, '\n');
    if (replace != NULL)
        replace[0] = '\0';


    return;
}
/* ------------------------------------------------------------------------- */






unsigned long CHARM(shc_read_gfc)(const char *pathname,
                                  unsigned long nmax,
                                  const char *epoch,
                                  CHARM(shc) *shcs,
                                  CHARM(err) *err)
{
    /* Open "pathname" to read */
    /* --------------------------------------------------------------------- */
    FILE *fptr = fopen(pathname, "r");
    char err_msg[CHARM_ERR_MAX_MSG];
    if (fptr == NULL)
    {
        sprintf(err_msg, "Couldn't open \"%s\".", pathname);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       err_msg);
        return CHARM_SHC_NMAX_ERROR;
    }
    /* --------------------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
    char line[SHC_READ_GFC_NLINE];
    char key_str[SHC_READ_GFC_NSTR];
    char val_str[SHC_READ_GFC_NSTR];
    REAL mu_file = PREC(0.0);
    REAL r_file  = PREC(0.0);
    unsigned long nmax_file = CHARM_SHC_NMAX_ERROR;
    /* --------------------------------------------------------------------- */






    /* Read the header of the "gfc" file */
    /* --------------------------------------------------------------------- */
    /* These variables are used to mark that a specific keyword of the "gfc"
     * file has been found */
    _Bool nmax_found, mu_found, r_found, errors_found, eoh_found, boh_found;
    nmax_found = mu_found = r_found = errors_found = eoh_found = boh_found = 0;


    /* These variables are set to true if the "begin_of_head" keyword was not
     * found before reaching the keyword and an error was encountered when
     * processing the value of that keyword.  This needed is because
     * "begin_of_head" is not a mandatory keyword of "gfc" files, in which case
     * the keywords have to be taken from the comments section.  But if
     * "begin_of_head" will appear later, the error in conversion is not a real
     * error, because it happened in the comments section. */
    _Bool nmax_err, mu_err, r_err, errors_err, norm_err, format_err;
    nmax_err = mu_err = r_err = errors_err = norm_err = format_err = 0;


    /* These variables are used to store the value of the keyword that couldn't
     * be processed if "nmax_err", "mu_err", ... are true */
    char nmax_err_str[SHC_READ_GFC_NSTR];
    char mu_err_str[SHC_READ_GFC_NSTR];
    char r_err_str[SHC_READ_GFC_NSTR];
    char errors_err_str[SHC_READ_GFC_NSTR];
    char norm_err_str[SHC_READ_GFC_NSTR];
    char format_err_str[SHC_READ_GFC_NSTR];


    int ret;


    /* Type of errors of "gfc" files. */
    _Bool errors_no, errors_cal, errors_form, errors_cal_form;
    errors_no = errors_cal = errors_form = errors_cal_form = 0;


    _Bool icgem1d0, icgem2d0;
    icgem1d0 = icgem2d0 = 0;


    const char *wrong_value_of_keyword    = NULL;
    const char *wrong_value_required_type = NULL;


    do
    {
        /* Read a line of the file */
        if (fgets(line, SHC_READ_GFC_NLINE, fptr) == NULL)
        {
            sprintf(err_msg, "Couldn't find the \"end_of_head\" keyword in "
                             "\"%s\" or couldn't read the file.", pathname);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                           err_msg);
            goto EXIT;
        }


        /* Get the first two entries of "line" to test for the keywords of
         * "gfc" files */
        errno = 0;
        ret = sscanf(line, "%s %s", key_str, val_str);
        if (errno)
        {
            sprintf(err_msg, "Couldn't read \"%s\" with \"sscanf\".", pathname);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                           err_msg);
            goto EXIT;
        }


        if ((ret == 0) || (ret == EOF))
            continue;


        /* Is this the line specifying the begin_of_head keyword? */
        if (strcmp(key_str, SHC_READ_GFC_BOH) == 0)
            boh_found = 1;
        else
            boh_found = 0;


        /* Is this the line specifying the end_of_head keyword? */
        if (strcmp(key_str, SHC_READ_GFC_EOH) == 0)
            eoh_found = 1;
        else
            eoh_found = 0;


        /* From now on, we are searching for a keyword and its value, so "ret"
         * must be two */
        if (ret != 2)
            continue;


        if (strcmp(key_str, SHC_READ_GFC_NMAX) == 0)
        {
            nmax_err  = 0;
            nmax_file = CHARM(misc_str2ul)(val_str, "", err);
            if (!CHARM(err_isempty)(err))
            {
                strcpy(nmax_err_str, val_str);
                nmax_err = 1;
                CHARM(err_reset)(err);
            }


            if (boh_found && nmax_err)
            {
FAILURE_NMAX:
                wrong_value_required_type = "unsigned long int";
                goto FAILURE_STR2NUM;
            }
            nmax_found = 1;
        }
        else if (strcmp(key_str, SHC_READ_GFC_GM) == 0)
        {
            mu_err   = 0;
            mu_file  = CHARM(misc_str2real)(val_str, "", err);
            if (!CHARM(err_isempty)(err))
            {
                strcpy(mu_err_str, val_str);
                mu_err = 1;
                CHARM(err_reset)(err);
            }


            if (boh_found && mu_err)
            {
FAILURE_MU:
                wrong_value_required_type = "REAL";
                goto FAILURE_STR2NUM;
            }
            mu_found = 1;
        }
        else if (strcmp(key_str, SHC_READ_GFC_R) == 0)
        {
            r_err   = 0;
            r_file  = CHARM(misc_str2real)(val_str, "", err);
            if (!CHARM(err_isempty)(err))
            {
                strcpy(r_err_str, val_str);
                r_err = 1;
                CHARM(err_reset)(err);
            }


            if (boh_found && r_err)
            {
FAILURE_R:
                wrong_value_required_type = "REAL";
                goto FAILURE_STR2NUM;
            }
            r_found = 1;
        }
        else if (strcmp(key_str, SHC_READ_GFC_ERRORS) == 0)
        {
            errors_err = 0;
            if (strcmp(val_str, SHC_READ_GFC_ERRORS_NO) == 0)
                errors_no = 1;
            else if (strcmp(val_str, SHC_READ_GFC_ERRORS_CALIBRATED) == 0)
                errors_cal = 1;
            else if (strcmp(val_str, SHC_READ_GFC_ERRORS_FORMAL) == 0)
                errors_form = 1;
            else if (strcmp(val_str,
                            SHC_READ_GFC_ERRORS_CALIBRATED_AND_FORMAL) == 0)
                errors_cal_form = 1;
            else
            {
                errors_err = 1;
                strcpy(errors_err_str, val_str);
                if (boh_found)
                {
FAILURE_ERRORS:
                    wrong_value_of_keyword = SHC_READ_GFC_ERRORS;
                    goto FAILURE_UNSUPPORTED_KEYWORD;
                }
            }
            errors_found = 1;
        }
        else if (strcmp(key_str, SHC_READ_GFC_NORM) == 0)
        {
            norm_err = 0;
            if (strcmp(val_str, SHC_READ_GFC_NORM_FULL) != 0)
            {
                norm_err = 1;
                strcpy(norm_err_str, val_str);
                if (boh_found)
                {
FAILURE_NORM:
                    wrong_value_of_keyword = SHC_READ_GFC_NORM;
                    goto FAILURE_UNSUPPORTED_KEYWORD;
                }
            }
            /* "norm" is an optional keyword, so no "norm_found" here */
        }
        else if (strcmp(key_str, SHC_READ_GFC_FORMAT) == 0)
        {
            format_err = 0;
            if (strcmp(val_str, SHC_READ_GFC_FORMAT1d0) == 0)
                icgem1d0 = 1;
            else if (strcmp(val_str, SHC_READ_GFC_FORMAT2d0) == 0)
                icgem2d0 = 1;
            else
            {
                format_err = 1;
                strcpy(format_err_str, val_str);
                if (boh_found)
                {
FAILURE_FORMAT:
                    wrong_value_of_keyword = SHC_READ_GFC_FORMAT;
                    goto FAILURE_UNSUPPORTED_KEYWORD;
                }
            }
            /* "format" is an optional keyword, so no "format_found" here */
        }

    } while (!eoh_found);


    /* Check whether all mandatory keywords were found */
    const char *not_found_keyword = NULL;
    if (!nmax_found)
        not_found_keyword = SHC_READ_GFC_NMAX;
    else if (!mu_found)
        not_found_keyword = SHC_READ_GFC_GM;
    else if (!r_found)
        not_found_keyword = SHC_READ_GFC_R;
    else if (!errors_found)
        not_found_keyword = SHC_READ_GFC_ERRORS;
    if (not_found_keyword != NULL)
    {
        sprintf(err_msg, "Couldn't find the \"%s\" keyword before reaching "
                         "the \"%s\" keyword in \"%s\".",
                not_found_keyword, SHC_READ_GFC_EOH, pathname);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                       CHARM_EFILEIO, err_msg);
        goto EXIT;
    }


    /* Now, after reading the comment and header sections, we have go back
     * a bit and to check whether the values of the individual keywords were
     * properly processed *if the "begin_of_head" keyword was *not* found*.
     * This is because if "begin_of_head" was not found, we had to take the
     * values from the comments section.  But in this case, we did not want to
     * throw errors if the values were invalid, simply because the
     * "begin_of_head" keyword might appear later, implying this was a comment
     * section indeed.  This complicated procedure results from the definition
     * of "gfc" files. */
    if (!boh_found)
    {
        if (nmax_err)
        {
            strcpy(val_str, nmax_err_str);
            strcpy(key_str, SHC_READ_GFC_NMAX);
            goto FAILURE_NMAX;
        }
        else if (mu_err)
        {
            strcpy(val_str, mu_err_str);
            strcpy(key_str, SHC_READ_GFC_GM);
            goto FAILURE_MU;
        }
        else if (r_err)
        {
            strcpy(val_str, r_err_str);
            strcpy(key_str, SHC_READ_GFC_R);
            goto FAILURE_R;
        }
        else if (errors_err)
        {
            strcpy(val_str, errors_err_str);
            strcpy(key_str, SHC_READ_GFC_ERRORS);
            goto FAILURE_ERRORS;
        }
        else if (norm_err)
        {
            strcpy(val_str, norm_err_str);
            strcpy(key_str, SHC_READ_GFC_NORM);
            goto FAILURE_NORM;
        }
        else if (format_err)
        {
            strcpy(val_str, format_err_str);
            strcpy(key_str, SHC_READ_GFC_FORMAT);
            goto FAILURE_FORMAT;
        }
    }


    if (CHARM(shc_read_nmax_only)(nmax, shcs))
        goto EXIT;


    shcs->mu = mu_file;
    shcs->r  = r_file;


    /* If "format" is not specified, default format is "icgem1.0" as per
     * definition of "gfc". */
    if (!icgem1d0 && !icgem2d0)
        icgem1d0 = 1;
    /* --------------------------------------------------------------------- */






    /* Check "epoch" and convert it to a "REAL_EPOCH" data type */
    /* --------------------------------------------------------------------- */
    /* "epoch = NULL" is a valid input to "epoch_fraction".  This case will be
     * treated later. */
    REAL_EPOCH t = epoch_fraction(epoch, err);
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
        sprintf(err_msg, "\"shcs->nmax = %lu\" cannot be smaller than the "
                         "input parameter \"nmax = %lu\".", shcs->nmax, nmax);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto EXIT;
    }


    if (nmax_file < nmax)
    {
        sprintf(err_msg, "Couldn't read coefficients up to degree "
                         "\"nmax = %lu\", because the maximum degree in "
                         "\"%s\" is \"%lu\" only.", nmax, pathname, nmax_file);
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       err_msg);
        goto EXIT;
    }


    if (icgem2d0 && (epoch == NULL))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"epoch\" cannot be \"NULL\" for \"gfc\" files of the "
                       "\"icgem2.0\" format.");
        goto EXIT;
    }
    /* --------------------------------------------------------------------- */






    /* Read the table of spherical harmonic coefficients */
    /* --------------------------------------------------------------------- */
    int ns;
    char s0[SHC_READ_GFC_NSTR];
    char s1[SHC_READ_GFC_NSTR];
    char s2[SHC_READ_GFC_NSTR];
    char s3[SHC_READ_GFC_NSTR];
    char s4[SHC_READ_GFC_NSTR];
    char s5[SHC_READ_GFC_NSTR];
    char s6[SHC_READ_GFC_NSTR];
    char s7[SHC_READ_GFC_NSTR];
    char s8[SHC_READ_GFC_NSTR];
    char s9[SHC_READ_GFC_NSTR];
    char s10[SHC_READ_GFC_NSTR];
    char s11[SHC_READ_GFC_NSTR];


    /* The "t0" and "t1" epochs loaded from the file as strings */
    char *t0_str = NULL;
    char *t1_str = NULL;


    /* The "t0" and "t1" epochs as floating points */
    REAL_EPOCH t0 = PREC(0.0);
    REAL_EPOCH t1 = PREC(0.0);


    /* The phase factor of the time variable coefficients as a floating point
     * */
    REAL p;


    /* The phase factor of the time variable coefficients as a string */
    char *p_str = NULL;


    /* Useful substitutions */
    const REAL twopi = PREC(2.0) * PI;
    REAL tmp;
    _Bool dot_found, trnd_found, asin_found, acos_found;
    dot_found = trnd_found = asin_found = acos_found = 0;


    /* If the format of the "gfc" file is "icgem1.0" and no epoch was specified
     * by the user, then "icgem1d0_no_epoch" is set to true.  Otherwise, it is
     * false.
     *
     * If "icgem1d0_no_epoch" is true, the default epoch for each coefficient
     * is used */
    _Bool icgem1d0_no_epoch = 0;
    if (icgem1d0 && (epoch == NULL))
        icgem1d0_no_epoch = 1;


    /* At first, reset all coefficients in "shcs" to zero. */
    CHARM(shc_reset_coeffs)(shcs);


    unsigned long n, m, n_tmp, m_tmp;
    n = m = n_tmp = m_tmp = 0;
    REAL cnm, snm;
    cnm = snm = PREC(0.0);
    while (fgets(line, SHC_READ_GFC_NLINE, fptr) != NULL)
    {
        /* Read line from the data section */
        errno = 0;
        ns = sscanf(line, "%s %s %s %s %s %s %s %s %s %s %s %s",
                    s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11);
        if (errno)
        {
            sprintf(err_msg, "Couldn't read with \"sscanf\" from \"%s\".",
                    pathname);
            CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                           CHARM_EFILEIO, err_msg);
            goto EXIT;
        }


        if (ns == EOF)
            /* Probably an empty line which is valid in "gfc" files */
            continue;


        /* Check for the keyword and process "line" accordingly */
        if (strcmp(s0, SHC_READ_GFC_GFC) == 0)
        {
            /* For "SHC_READ_GFC_GFC", the number of allowed values is 4, 6 or
             * 8 according to the specification (5, 7 or 9 after adding the gfc
             * keyword itself).  Interestingly, some older models available
             * from "ICGEM" with "errors_no == 1" may have even only three
             * values per line if "m == 0" which is why we accept "ns" to be
             * also "4".  At this point, we do not know whether "m == 0" or "m
             * != 0", so we have a wildcard for "ns != 4".  We do the check
             * later because the value of "m" in "line" has not yet been
             * determined. */
            if ((errors_no && (ns < 4)) ||
                ((errors_cal || errors_form) && (ns < 7)) ||
                (errors_cal_form && (ns < 9)))
            {
                WRONG_NUMBER_OF_ENTRIES;
            }


            READ_DEG_ORD(n, s1, "degree", pathname);
            if (n > nmax)
                continue;
            READ_DEG_ORD(m, s2, "order", pathname);
            READ_CNM_SNM(cnm, s3, SHC_READ_GFC_GFC, n, m, pathname);
            if (m == 0)
                snm = PREC(0.0);
            else
            {
                if (ns < 5)
                {
                    WRONG_NUMBER_OF_ENTRIES;
                }


                READ_CNM_SNM(snm, s4, SHC_READ_GFC_GFC, n, m, pathname);
            }
        }
        else if (strcmp(s0, SHC_READ_GFC_GFCT) == 0)
        {
            if (icgem1d0)
            {
                if ((errors_no && (ns < 6)) ||
                    ((errors_cal || errors_form) && (ns < 8)) ||
                    (errors_cal_form && (ns < 10)))
                {
                    WRONG_NUMBER_OF_ENTRIES;
                }
            }
            else if (icgem2d0)
            {
                if ((errors_no && (ns < 7)) ||
                    ((errors_cal || errors_form) && (ns < 9)) ||
                    (errors_cal_form && (ns < 11)))
                {
                    WRONG_NUMBER_OF_ENTRIES;
                }
            }


            /* Get the degree and order values */
            /* ------------------------------------------------------------- */
            READ_DEG_ORD(n, s1, "degree", pathname);
            if (n > nmax)
                continue;
            READ_DEG_ORD(m, s2, "order", pathname);
            /* ------------------------------------------------------------- */


            /* Get the epochs */
            /* ------------------------------------------------------------- */
            if (errors_no)
                t0_str = s5;
            else if (errors_cal || errors_form)
                t0_str = s7;
            else if (errors_cal_form)
                t0_str = s9;


            t0 = epoch_fraction(t0_str, err);
            if (!CHARM(err_isempty)(err))
            {
                CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                goto EXIT;
            }


            if (icgem2d0)
            {
                if (errors_no)
                    t1_str = s6;
                else if (errors_cal || errors_form)
                    t1_str = s8;
                else if (errors_cal_form)
                    t1_str = s10;


                t1 = epoch_fraction(t1_str, err);
                if (!CHARM(err_isempty)(err))
                {
                    CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                    goto EXIT;
                }


                if (!is_within_time_period(t, t0, t1))
                    continue;
            }
            /* ------------------------------------------------------------- */


            /* This code block needs to be placed after the check
             * "is_within_time_period" */
            if (trnd_found && dot_found)
            {
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFILEIO, "");
                goto EXIT;
            }
            else if ((trnd_found || dot_found) && asin_found && acos_found)
                /* All keywords associated with "gfct" were found, so reset
                 * reset some variables to false, so that we can continue
                 * reading */
                trnd_found = dot_found = asin_found = acos_found = 0;
            else if (trnd_found || dot_found || asin_found || acos_found)
            {
                remove_new_line_char(line);
                sprintf(err_msg, "Wrong format of \"%s\".  At least one of "
                        "the following keywords was found before reaching the "
                        "\"%s\" keyword or is missing: \"%s\", \"%s\", "
                        "\"%s\", \"%s\".  Stopped reading at line: \"%s\".",
                        pathname, SHC_READ_GFC_GFCT, SHC_READ_GFC_TRND,
                        SHC_READ_GFC_DOT, SHC_READ_GFC_ASIN,
                        SHC_READ_GFC_ACOS, line);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFILEIO, err_msg);
                goto EXIT;
            }


            /* Get the coefficients */
            READ_CNM_SNM(cnm, s3, SHC_READ_GFC_GFCT, n, m, pathname);
            READ_CNM_SNM(snm, s4, SHC_READ_GFC_GFCT, n, m, pathname);


            trnd_found = dot_found = asin_found = acos_found = 0;
        }
        else if ((strcmp(s0, SHC_READ_GFC_TRND) == 0) ||
                 (strcmp(s0, SHC_READ_GFC_DOT) == 0))
        {
            if (icgem1d0)
            {
                if ((errors_no && (ns < 5)) ||
                    ((errors_cal || errors_form) && (ns < 7)) ||
                    (errors_cal_form && (ns < 9)))
                {
                    WRONG_NUMBER_OF_ENTRIES;
                }
            }
            else if (icgem2d0)
            {
                if ((errors_no && (ns < 7)) ||
                    ((errors_cal || errors_form) && (ns < 9)) ||
                    (errors_cal_form && (ns < 11)))
                {
                    WRONG_NUMBER_OF_ENTRIES;
                }
            }


            /* Get the degree and order values */
            /* ------------------------------------------------------------- */
            READ_DEG_ORD(n_tmp, s1, "degree", pathname);
            if (n_tmp > nmax)
                continue;
            READ_DEG_ORD(m_tmp, s2, "order", pathname);
            EXPECTED_DEGREE_ORDER(s0, n, m, pathname, n_tmp, m_tmp, s3, s4);
            /* ------------------------------------------------------------- */


            /* Get the epoch.  For the "icgem1.0" format, the "t0" value is
             * taken from the previous loop run.  For "icgem2.0", the "t0" and
             * "t1" values are taken from the current line of the file. */
            if (icgem2d0)
            {
                if (errors_no)
                {
                    t0_str = s5;
                    t1_str = s6;
                }
                else if (errors_cal || errors_form)
                {
                    t0_str = s7;
                    t1_str = s8;
                }
                else if (errors_cal_form)
                {
                    t0_str = s9;
                    t1_str = s10;
                }


                t0 = epoch_fraction(t0_str, err);
                if (!CHARM(err_isempty)(err))
                {
                    CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                    goto EXIT;
                }
                t1 = epoch_fraction(t1_str, err);
                if (!CHARM(err_isempty)(err))
                {
                    CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                    goto EXIT;
                }


                if (!is_within_time_period(t, t0, t1))
                    continue;
            }


            if (strcmp(s0, SHC_READ_GFC_DOT) == 0)
            {
                READ_CNM_SNM(cnm, s3, SHC_READ_GFC_DOT, n, m, pathname);
                READ_CNM_SNM(snm, s4, SHC_READ_GFC_DOT, n, m, pathname);
            }
            else
            {
                READ_CNM_SNM(cnm, s3, SHC_READ_GFC_TRND, n, m, pathname);
                READ_CNM_SNM(snm, s4, SHC_READ_GFC_TRND, n, m, pathname);
            }


            if (icgem1d0_no_epoch)
                t = t0;
            tmp  = epoch_diff(t, t0);
            cnm *= tmp;
            snm *= tmp;


            if (strcmp(s0, SHC_READ_GFC_DOT) == 0)
            {
                /* The "dot" keyword is not associated with the "asin" and
                 * "acos" terms.  Here, this can be easily solved if we pretend
                 * that we have already found "asin" and "acos", so that we can
                 * continue in reading the file */
                dot_found = asin_found = acos_found = 1;
            }
            else
                trnd_found = 1;
        }
        else if ((strcmp(s0, SHC_READ_GFC_ASIN) == 0) ||
                 (strcmp(s0, SHC_READ_GFC_ACOS) == 0))
        {
            if (icgem1d0)
            {
                if ((errors_no && (ns < 6)) ||
                    ((errors_cal || errors_form) && (ns < 8)) ||
                    (errors_cal_form && (ns < 10)))
                {
                    WRONG_NUMBER_OF_ENTRIES;
                }
            }
            else if (icgem2d0)
            {
                if ((errors_no && (ns < 8)) ||
                    ((errors_cal || errors_form) && (ns < 10)) ||
                    (errors_cal_form && (ns < 12)))
                {
                    WRONG_NUMBER_OF_ENTRIES;
                }
            }


            /* Get the degree and order values */
            /* ------------------------------------------------------------- */
            READ_DEG_ORD(n_tmp, s1, "degree", pathname);
            if (n_tmp > nmax)
                continue;
            READ_DEG_ORD(m_tmp, s2, "order", pathname);
            EXPECTED_DEGREE_ORDER(s0, n, m, pathname, n_tmp, m_tmp, s3, s4);
            /* ------------------------------------------------------------- */


            /* Get the epoch.  For the "icgem1.0" format, the "t0" value is
             * taken from the previous loop run.  For "icgem2.0", the "t0" and
             * "t1" values are taken from the current line of the file. */
            if (icgem2d0)
            {
                if (errors_no)
                {
                    t0_str = s5;
                    t1_str = s6;
                }
                else if (errors_cal || errors_form)
                {
                    t0_str = s7;
                    t1_str = s8;
                }
                else if (errors_cal_form)
                {
                    t0_str = s9;
                    t1_str = s10;
                }


                t0 = epoch_fraction(t0_str, err);
                if (!CHARM(err_isempty)(err))
                {
                    CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                    goto EXIT;
                }
                t1 = epoch_fraction(t1_str, err);
                if (!CHARM(err_isempty)(err))
                {
                    CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
                    goto EXIT;
                }


                if (!is_within_time_period(t, t0, t1))
                    continue;
            }


            if (strcmp(s0, SHC_READ_GFC_ASIN) == 0)
            {
                READ_CNM_SNM(cnm, s3, SHC_READ_GFC_ASIN, n, m, pathname);
                READ_CNM_SNM(snm, s4, SHC_READ_GFC_ASIN, n, m, pathname);
            }
            else
            {
                READ_CNM_SNM(cnm, s3, SHC_READ_GFC_ACOS, n, m, pathname);
                READ_CNM_SNM(snm, s4, SHC_READ_GFC_ACOS, n, m, pathname);
            }


            /* Get the phase factor */
            if (icgem1d0)
            {
                if (errors_no)
                    p_str = s5;
                else if (errors_cal || errors_form)
                    p_str = s7;
                else if (errors_cal_form)
                    p_str = s9;
            }
            else if (icgem2d0)
            {
                if (errors_no)
                    p_str = s7;
                else if (errors_cal || errors_form)
                    p_str = s9;
                else if (errors_cal_form)
                    p_str = s11;
            }
            p = CHARM(misc_str2real)(p_str, "", err);
            if (!CHARM(err_isempty)(err))
            {
                sprintf(err_msg, "Failed to convert the phase factor \"%s\" "
                                 "of time variable coefficients of degree "
                                 "\"%lu\" and order \"%lu\" from \"%s\" "
                                 "to a \"REAL\" data format.", p_str, n, m,
                                 pathname);
                CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                               CHARM_EFILEIO, err_msg);
                goto EXIT;
            }


            if (icgem1d0_no_epoch)
                t = t0;
            tmp  = epoch_diff(t, t0);
            tmp = (twopi / p) * tmp;
            if (strcmp(s0, SHC_READ_GFC_ASIN) == 0)
            {
                tmp        = SIN(tmp);
                asin_found = 1;
            }
            else
            {
                tmp        = COS(tmp);
                acos_found = 1;
            }
            cnm *= tmp;
            snm *= tmp;
        }
        else
            /* Lines starting with any other keyword are comments, so should be
             * skipped */
            continue;


        shcs->c[m][n - m] += cnm;
        shcs->s[m][n - m] += snm;
    }
    /* --------------------------------------------------------------------- */


EXIT:
    fclose(fptr);
    return nmax_file;


FAILURE_STR2NUM:
    CHARM(err_reset)(err);
    sprintf(err_msg, "Failed to convert the \"%s\" value of the \"%s\" "
                     "keyword from \"%s\" to the \"%s\" "
                     "data format.", val_str, key_str, pathname,
                     wrong_value_required_type);
    CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO, err_msg);
    goto EXIT;


FAILURE_UNSUPPORTED_KEYWORD:
    sprintf(err_msg, "\"%s\" is not a supported value of the \"%s\" keyword "
                     "in \"%s\".", val_str, wrong_value_of_keyword, pathname);
    CHARM(err_set)(err, __FILE__, __LINE__, __func__,
                   CHARM_EFILEIO, err_msg);
    goto EXIT;
}
