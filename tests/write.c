/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "write_array.h"
#include "write.h"
/* ------------------------------------------------------------------------- */






/* Writes "n" data elements from the array pointed to by "f" to a file pointed
 * to by "file".  */
/* ------------------------------------------------------------------------- */
long int write(char *file, REAL *f, size_t n)
{
    /* Read the reference data from "file" */
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
    long int e = write_array(f, n, fptr);


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
