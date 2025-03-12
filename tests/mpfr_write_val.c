/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "mpfr_write_val.h"
/* ------------------------------------------------------------------------- */






/* Write "mpfr_t" to a stream.  On success, returned is "0".  Otherwise, "1" is
 * returned. */
/* ------------------------------------------------------------------------- */
int mpfr_write_val(mpfr_t val,
                   FILE *fptr)
{
    size_t ret = mpfr_out_str(fptr, 10, 0, val, MPFR_RNDN);
    fprintf(fptr, "\n");


    return (ret > 0) ? 0 : 1;
}
/* ------------------------------------------------------------------------- */

