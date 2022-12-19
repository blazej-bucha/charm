/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(misc_str2real)(char str[], const char *err_msg, CHARM(err) *err)
{
    char *end_ptr;
    errno = 0;


    /* The "strtod", etc., functions that convert strings to floating point
     * numbers do not consider the Fortran's "d" and "D" decimal exponents.
     * Here, we therefore check for the presence of the first occurrence of the
     * 'D' and 'd' characters in "str".  If there is a match, the character
     * found is replaced by 'E' or 'e', respectively.  It is sufficient to
     * check for the first presence of the characters, as more than one 'D' or
     * 'd' implies as incorrect format anyway. */
    char *match = strchr(str, 'D');
    if (match != NULL)
        *match = 'E';
    match = strchr(str, 'd');
    if (match != NULL)
        *match = 'e';


    REAL r = STR2REAL(str, &end_ptr);


    if ((end_ptr == str) || (errno != 0))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       err_msg);


    return r;
}
