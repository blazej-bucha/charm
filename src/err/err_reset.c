/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <string.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(err_reset)(CHARM(err) *err)
{
    if (err == NULL)
        return;


    err->level = 0;


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        memset(err->file[e], '\0', CHARM_ERR_MAX_FILE * sizeof(char));


    memset(err->line, 0, CHARM_ERR_MAX_LEVEL * sizeof(size_t));


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        memset(err->func[e], '\0', CHARM_ERR_MAX_FUNC * sizeof(char));


    err->code = CHARM_SUCCESS;


    memset(err->msg, '\0', CHARM_ERR_MAX_MSG * sizeof(char));


    err->saturated = 0;


    return;
}
