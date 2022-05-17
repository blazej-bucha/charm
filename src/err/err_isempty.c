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


    if (err->level != 0)
        return 0;


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        for (size_t i = 0; i < CHARM_ERR_MAX_FILE; i++)
            if (err->file[e][i] != '\0')
                return 0;


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        if (err->line[e] != 0)
            return 0;


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        for (size_t i = 0; i < CHARM_ERR_MAX_FUNC; i++)
            if (err->func[e][i] != '\0')
                return 0;


    if (err->code != CHARM_SUCCESS)
        return 0;


    for (size_t i = 0; i < CHARM_ERR_MAX_MSG; i++)
        if (err->msg[i] != '\0')
            return 0;


    return 1;
}
