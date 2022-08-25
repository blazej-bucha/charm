/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "../prec.h"
#include "../err/err_set.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(misc_str2real)(const char *str, const char *err_msg,
                          CHARM(err) *err)
{
    char *end_ptr;
    errno = 0;
    REAL r = STR2REAL(str, &end_ptr);


    if ((end_ptr == str) || (errno != 0))
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFILEIO,
                       err_msg);


    return r;
}
