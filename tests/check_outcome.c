/* Header files */
/* ------------------------------------------------------------------------- */
#include <stdio.h>
/* ------------------------------------------------------------------------- */






void check_outcome(long int err)
{
    if (err)
        printf("\n            %ld possibly wrong result%s\n",
               err, (err == 1) ? "" : "s");
    else
        printf("ok\n");
}
