/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "../src/misc/misc_fprintf_real.h"
#include "write_val.h"
/* ------------------------------------------------------------------------- */






/* Write floating point value to a stream.  On success, returned is "0".
 * Otherwise, "1" is returned. */
/* ------------------------------------------------------------------------- */
int write_val(REAL val, FILE *fptr)
{
    long int ret = CHARM(misc_fprintf_real)(fptr, REAL_PRINT_FORMAT, val);
    fprintf(fptr, "%s", "\n");


    return (ret > 0) ? 0 : 1;
}
/* ------------------------------------------------------------------------- */

