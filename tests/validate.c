/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "cmp_arrays.h"
#include "../src/prec.h"
#include "../src/misc/misc_str2real.h"
/* ------------------------------------------------------------------------- */





#define NCHARS 256






/* Compares "n" elements from an array pointed to be "f" with respect to
 * reference data from "file" using the "eps" threshold. */
/* ------------------------------------------------------------------------- */
int validate(char *file, REAL *f, size_t n, REAL eps)
{
    /* Read the reference data from "file" */
    /* --------------------------------------------------------------------- */
    FILE *fptr = fopen(file, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "Failed to open \"%s\"\n", file);
        exit(CHARM_FAILURE);
    }


    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
    {
        fprintf(stderr, "Failed to initialize the \"err\" structure.\n");
        exit(CHARM_FAILURE);
    }


    /* A string to be loaded from the input file */
    char str[NCHARS];
    int num_entries;
    REAL entry_d;


    REAL *fref = (REAL *)malloc(n * sizeof(REAL));
    if (fref == NULL)
    {
        fprintf(stderr, "malloc failure.\n");
        exit(CHARM_FAILURE);
    }


    for (size_t i = 0; i < n; i++)
    {
        num_entries = fscanf(fptr, "%s", str);
        if (num_entries == EOF)
        {
            fprintf(stderr, "Too few data in \"%s\".", file);
            exit(1);
        }
        if (num_entries < 1)
        {
            fprintf(stderr, "Failed to read an entry from \"%s\".", file);
            exit(1);
        }


        entry_d = CHARM(misc_str2real)(str, "", err);
        if (!CHARM(err_isempty)(err))
        {
            fprintf(stderr, "Failed to convert \"%s\" to a floating point "
                            "data type.", str);
            exit(1);
        }


        fref[i] = entry_d;
    }


    fclose(fptr);
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */






    /* Compare "f" with respect to "fref" */
    /* --------------------------------------------------------------------- */
    int errnum = cmp_arrays(f, fref, n, eps);
    free(fref);


    if (errnum)
        printf("        Reference file \"%s\"\n", file);
    /* --------------------------------------------------------------------- */






    return errnum;
}
/* ------------------------------------------------------------------------- */