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
                          unsigned int line,
                          const char *func)
{
    if ((err == NULL) || err->saturated)
        return;


    /* Get the error level to which we will write the data. */
    unsigned int curr_level = err->level;


    /* Write the file name.  If the length of the string in "file" is longer
     * than "CHARM_ERR_MAX_FILE" characters, some information will be lost.
     * But at least a buffer overflow should never happen.  At any rate,
     * "CHARM_ERR_MAX_FILE" is large enough, so we should really be on the safe
     * side in any reasonable circumstances.  The same applies to other parts
     * of this function. */
    strncpy(err->file[curr_level], file, CHARM_ERR_MAX_FILE - 1);
    err->file[curr_level][CHARM_ERR_MAX_FILE - 1] = '\0';


    /* Write the line number */
    err->line[curr_level] = line;


    /* Write the function name */
    strncpy(err->func[curr_level], func, CHARM_ERR_MAX_FUNC - 1);
    err->func[curr_level][CHARM_ERR_MAX_FUNC - 1] = '\0';


    /* Finally, increase "err->level" by "1", since we have just written some
     * stuff to a new error level. */
    CHARM(err_inc_level)(err);


    return;
}
