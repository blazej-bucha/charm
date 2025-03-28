/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <string.h>
#include "../prec.h"
#include "err_inc_level.h"
#include "err_set.h"
/* ------------------------------------------------------------------------- */






void CHARM(err_set)(CHARM(err) *err,
                    const char *file,
                    size_t line,
                    const char *func,
                    int code,
                    const char *msg)
{
    if ((err == NULL) || (err->saturated))
        return;


    /* Get the error level to which we will write the data. */
    size_t curr_level = err->level;


    /* Write the file name. */
    snprintf(err->file[curr_level], CHARM_ERR_MAX_FILE, file);


    /* Write the line number */
    err->line[curr_level] = line;


    /* Write the function name */
    snprintf(err->func[curr_level], CHARM_ERR_MAX_FUNC, func);


    /* Write the error code */
    err->code = code;


    /* Write the error message */
    snprintf(err->msg, CHARM_ERR_MAX_MSG, msg);


    /* Finally, increase "err->level" by "1", since we have just written some
     * stuff to a new error level. */
    CHARM(err_inc_level)(err);


    return;
}
