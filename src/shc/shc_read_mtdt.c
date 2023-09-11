/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../misc/misc_str2ul.h"
#include "../misc/misc_str2real.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
/* ------------------------------------------------------------------------- */






/* Symbolic constants */
/* ------------------------------------------------------------------------- */
#undef NSTR
#define NSTR (128)


#undef NLINE
#define NLINE (2048)
/* ------------------------------------------------------------------------- */






/* Reads metadata of a "shc" structure from a text file: the maximum harmonic
 * degree, the scaling constant and the radius of the reference sphere.  This
 * function is not a part of API. */
void CHARM(shc_read_mtdt)(FILE *stream,
                          unsigned long *nmax,
                          REAL *mu,
                          REAL *r,
                          CHARM(err) *err)
{
    /* Char arrays to store the strings to be loaded from the input file */
    /* --------------------------------------------------------------------- */
    char line[NLINE];
    char nmax_str[NSTR];
    char mu_str[NSTR];
    char r_str[NSTR];
    /* --------------------------------------------------------------------- */






    /* Get the first line of the input file */
    /* --------------------------------------------------------------------- */
    if (fgets(line, NLINE, stream) == NULL)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to read the first line of the input file, "
                       "where the maximum harmonic degree, the scaling "
                       "constant and the radius of the reference sphere "
                       "are supposed to be given.");
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Read the metadata of the "shcs" structure */
    /* --------------------------------------------------------------------- */
    /* Get the maximum harmonic degree, the scaling constant and the radius of
     * the reference sphere as a string */
    /* ..................................................................... */
    if (sscanf(line, "%s %s %s", nmax_str, mu_str, r_str) != 3)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "Failed to read the metadata from the input file "
                       "(the first line of the file).");
        return;
    }
    /* ..................................................................... */


    /* Convert the strings to numerical data formats */
    /* ..................................................................... */
    *nmax = CHARM(misc_str2ul)(nmax_str, "Failed to convert the maximum "
                               "harmonic degree to the "
                               "\"unsigned long int\" data format.", err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    *mu = CHARM(misc_str2real)(mu_str, "Failed to convert the scaling "
                               "parameter to the \"REAL\" data format.", err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    *r = CHARM(misc_str2real)(r_str, "Failed to convert the radius of the "
                              "reference sphere to the \"REAL\" data format.",
                              err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    if (*r <= PREC(0.0))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       "The radius of the reference "
                       "sphere must be larger than zero.");
        return;
    }
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */






    return;
}
