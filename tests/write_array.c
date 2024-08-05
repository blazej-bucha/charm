/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "write_val.h"
#include "write_array.h"
/* ------------------------------------------------------------------------- */






/* Write "n" floating point values from an array pointed to by "arr" to
 * a stream.  On success, returned is "0".  Any return value other than "0"
 * implies an error. */
/* ------------------------------------------------------------------------- */
long int write_array(REAL *arr, size_t n, FILE *fptr)
{
    long int ret = 0;


    for (size_t i = 0; i < n; i++)
        ret += write_val(arr[i], fptr);


    return ret;
}
/* ------------------------------------------------------------------------- */
