/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






_Bool CHARM(err_isempty)(const CHARM(err) *err)
{
    if (err == NULL)
        return 0;


    if (err->code != CHARM_SUCCESS)
        return 0;


    return 1;
}
