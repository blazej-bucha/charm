/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "parameters.h"
#include "error_messages.h"
#include "../src/mpfr/mpfr_flush_unreleased_memory.h"
#include "mpfr_cmp_arrays.h"
#include "mpfr_validate.h"
/* ------------------------------------------------------------------------- */






/* Compares "n" elements from an "mpfr_t" array pointed to by "f" with respect
 * to reference data from "file" using the "eps" threshold. */
/* ------------------------------------------------------------------------- */
long int mpfr_validate(char *file,
                       mpfr_t *f,
                       size_t n,
                       mpfr_t eps)
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
        fprintf(stderr, ERR_MSG_ERR);
        exit(CHARM_FAILURE);
    }


    /* A string to be loaded from the input file */
    char str[NSTR_LONG];
    int num_entries;


    mpfr_t *fref = (mpfr_t *)malloc(n * sizeof(mpfr_t));
    if (fref == NULL)
    {
        fprintf(stderr, CHARM_ERR_MALLOC_FAILURE"\n");
        exit(CHARM_FAILURE);
    }
    for (size_t i = 0; i < n; i++)
        mpfr_init2(fref[i], GFM_CAP_NBITS);


    for (size_t i = 0; i < n; i++)
    {
        num_entries = fscanf(fptr, "%s", str);
        if (num_entries == EOF)
        {
            fprintf(stderr, "Too few data in \"%s\".", file);
            exit(CHARM_FAILURE);
        }
        else if (num_entries < 1)
        {
            fprintf(stderr, "Failed to read an entry from \"%s\".", file);
            exit(CHARM_FAILURE);
        }


        mpfr_set_str(fref[i], str, 10, MPFR_RNDN);
    }


    fclose(fptr);
    CHARM(err_free)(err);
    /* --------------------------------------------------------------------- */


    /* Compare "f" with respect to "fref" */
    /* --------------------------------------------------------------------- */
    long int e = mpfr_cmp_arrays(f, fref, n, eps);
    if (e)
        printf("\n        Reference file \"%s\"\n", file);


    for (size_t i = 0; i < n; i++)
        mpfr_clear(fref[i]);
    free(fref);
    /* --------------------------------------------------------------------- */


    FLUSH_UNRELEASED_MEMORY;


    return e;
}
/* ------------------------------------------------------------------------- */
