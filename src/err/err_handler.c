/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(err_handler)(CHARM(err) *err, _Bool terminate)
{
    if (CHARM(err_isempty)(err) || (err == NULL))
        return;


    /* Print details on the error to "stderr". */
    /* --------------------------------------------------------------------- */
    fprintf(stderr, "\n"
                    "-----------\n"
                    "CHarm ERROR\n"
                    "-----------\n"
                    "Error code: %u                   Traceback (most "
                    "recent call last)\n\n", err->code);


    for (size_t e = err->level - 1; e != (size_t)(-1); e--)
        fprintf(stderr, "    File \"%s\", line: %zu, function: \"%s\"\n\n",
                        err->file[e], err->line[e], err->func[e]);


    if (err->saturated)
        fprintf(stderr, "    Warning: The error structure is "
                        "saturated.  Most recent function calls may therefore "
                        "not be reported.\n\n");


    fprintf(stderr, "Error message: %s\n", err->msg);
    fprintf(stderr, "-----------\n");
    /* --------------------------------------------------------------------- */


    /* Terminate the program if the user decided so */
    if (terminate)
        exit(CHARM_FAILURE);


    /* And now reset "err" to the default values, so that it can be reused
     * again. */
    CHARM(err_reset)(err);
}
