/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>
#include "../src/prec.h"
#include "mpfr_write_array.h"
#include "mpfr_array2file.h"
/* ------------------------------------------------------------------------- */






/* Writes "n" data elements from the array pointed to by "f" to a file pointed
 * to by "file".  */
/* ------------------------------------------------------------------------- */
long int mpfr_array2file(char *file,
                         mpfr_t *f,
                         size_t n)
{
    /* Open "file" */
    /* --------------------------------------------------------------------- */
    FILE *fptr = fopen(file, "w");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to create \"%s\"\n", file);
        exit(CHARM_FAILURE);
    }
    /* --------------------------------------------------------------------- */


    /* Write "f" to "fptr" */
    /* --------------------------------------------------------------------- */
    long int e = mpfr_write_array(f, n, fptr);
    if (e)
    {
        fprintf(stderr, "Failed to write %lu reference values.\n", e);
        exit(CHARM_FAILURE);
    }


    fclose(fptr);
    /* --------------------------------------------------------------------- */


    return e;
}
/* ------------------------------------------------------------------------- */
