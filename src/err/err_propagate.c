/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../prec.h"
#include "err_inc_level.h"
#include "err_propagate.h"
/* ------------------------------------------------------------------------- */






void CHARM(err_propagate)(CHARM(err) *err,
                          const char *file,
                          size_t line,
                          const char *func)
{
    if ((err == NULL) || err->saturated)
        return;


    /* Get the error level to which we will write the data. */
    size_t curr_level = err->level;


    /* Write the file name. */
    snprintf(err->file[curr_level], CHARM_ERR_MAX_FILE, "%s", file);


    /* Write the line number */
    err->line[curr_level] = line;


    /* Write the function name */
    snprintf(err->func[curr_level], CHARM_ERR_MAX_FUNC, "%s", func);


    /* Finally, increase "err->level" by "1", since we have just written some
     * stuff to a new error level. */
    CHARM(err_inc_level)(err);


    return;
}
