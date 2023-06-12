/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/misc/misc_fprintf_real.h"
/* ------------------------------------------------------------------------- */






/* Write floating point value to a stream.  On success, returned is "0".
 * Otherwise, "1" is returned. */
/* ------------------------------------------------------------------------- */
int write_val(REAL val, FILE *fptr)
{
    char *format =
#if defined(CHARM_FLOAT)
                         "%0.9e";
#elif defined(CHARM_QUAD)
                         "%0.36Qe";
#else
                         "%0.18e";
#endif
    long int ret = CHARM(misc_fprintf_real)(fptr, format, val);
    fprintf(fptr, "\n");


    return (ret > 0) ? 0 : 1;
}
/* ------------------------------------------------------------------------- */

