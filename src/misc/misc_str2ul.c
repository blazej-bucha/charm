/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "misc_str2ul.h"
/* ------------------------------------------------------------------------- */






unsigned long CHARM(misc_str2ul)(const char *str, const char *err_msg,
                                 CHARM(err) *err)
{
    char *end_ptr;
    errno = 0;
    unsigned long ul = strtoul(str, &end_ptr, 10);


    /* "strtoul" silently converts negative values in "str" to "unsigned long",
     * so we also have to check "str" for the presence of the minus sign.  If
     * any, report an error. */
    if ((*end_ptr != '\0') || errno || (strchr(str, '-') != NULL))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       err_msg);


    return ul;
}
