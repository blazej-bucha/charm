/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "mpfr_write_val.h"
#include "mpfr_write_array.h"
/* ------------------------------------------------------------------------- */






/* Write "n" "mpfr_t" values from an array pointed to by "arr" to a stream.  On
 * success, returned is "0".  Any return value other than "0" implies an
 * error. */
/* ------------------------------------------------------------------------- */
long int mpfr_write_array(mpfr_t *arr,
                          size_t n,
                          FILE *fptr)
{
    long int ret = 0;


    for (size_t i = 0; i < n; i++)
        ret += mpfr_write_val(arr[i], fptr);


    return ret;
}
/* ------------------------------------------------------------------------- */
