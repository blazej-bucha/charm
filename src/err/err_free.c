/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(err_free)(CHARM(err) *err)
{
    if (err == NULL)
        return;


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        free(err->file[e]);


    for (size_t e = 0; e < CHARM_ERR_MAX_LEVEL; e++)
        free(err->func[e]);


    free(err->file);
    free(err->func);
    free(err->line);
    free(err->msg);
    free(err);


    return;
}
