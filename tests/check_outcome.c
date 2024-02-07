/* Header files */
/* ------------------------------------------------------------------------- */
#include <stdio.h>
#include "check_outcome.h"
/* ------------------------------------------------------------------------- */






void check_outcome(long int err)
{
    if (err)
        printf("\n            %ld possibly wrong result%s\n\n",
               err, (err == 1) ? "" : "s");
    else
        printf("ok\n");
}

