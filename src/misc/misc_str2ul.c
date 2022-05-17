/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






unsigned long CHARM(misc_str2ul)(char str[], char *err_msg, CHARM(err) *err)
{
    char *end_ptr;
    errno = 0;
    unsigned long ul = strtoul(str, &end_ptr, 10);


    if ((end_ptr == str) || (errno != 0))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       err_msg);


    return ul;
}
