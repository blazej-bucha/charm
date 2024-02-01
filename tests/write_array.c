/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "write_val.h"
#include "write_array.h"
/* ------------------------------------------------------------------------- */






/* Compares elements of two arrays up to some threshold */
/* ------------------------------------------------------------------------- */
long int write_array(REAL *arr, size_t n, FILE *fptr)
{
    long int ret = 0;


    for (size_t i = 0; i < n; i++)
        ret += write_val(arr[i], fptr);


    return ret;
}
/* ------------------------------------------------------------------------- */
