/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "misc_str2real.h"
/* ------------------------------------------------------------------------- */






/* Converts "str" to "REAL".  If failed, error message "err_msg" is entered to
 * "err".
 *
 * IMPORTANT NOTE: The function assumes write access to "str" if "d" or "D" is
 * used as decimal exponents.  "str" should thus not be a string literal,
 * etc. */
REAL CHARM(misc_str2real)(char *str,
                          const char *err_msg,
                          CHARM(err) *err)
{
    /* The "strtod", etc., functions that convert strings to floating point
     * numbers do not consider the Fortran's "d" and "D" decimal exponents.
     * Here, we therefore check for the presence of the first occurrence of the
     * 'D' and 'd' characters in "str".  If there is a match, the character
     * found is replaced by 'E' or 'e', respectively.  It is sufficient to
     * check for the first presence of the characters, as more than one 'D' or
     * 'd' implies an incorrect format anyway. */
    char *match_D = strchr(str, 'D');
    if (match_D != NULL)
        *match_D = 'E';
    char *match_d = strchr(str, 'd');
    if (match_d != NULL)
        *match_d = 'e';


    char *end_ptr;
    errno = 0;
    REAL r = STR2REAL(str, &end_ptr);


    if ((end_ptr == str) || errno)
        goto FAILURE;


    /* Any character in "str" after the last character used in the conversion
     * by "STR2REAL" is an error, except for one or more spaces or the
     * terminating null byte. */
    while (*end_ptr != '\0')
        if (!isspace(*end_ptr++))
            goto FAILURE;


EXIT:
    /* Before leaving the routine, set back "E" and "e" to "D" and "d" if
     * needed */
    if (match_D != NULL)
        *match_D = 'D';
    if (match_d != NULL)
        *match_d = 'd';


    return r;


FAILURE:
    CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO, err_msg);
    goto EXIT;
}
