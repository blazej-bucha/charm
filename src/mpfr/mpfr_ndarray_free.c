/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include <stdio.h>
#include <stdlib.h>
#include "mpfr_ndarray.h"
#include "mpfr_flush_unreleased_memory.h"
/* ------------------------------------------------------------------------- */






/* Internal function to free an "mpfr_ndarray" "x". */
void CHARM(mpfr_ndarray_free)(mpfr_ndarray *x)
{
    if (x == NULL)
        return;


    if ((x->owner) && (x->data != NULL))
        for (size_t i = 0; i < x->size; i++)
            mpfr_clear(x->data[i]);


    free(x->data);
    free(x->shape);
    free(x);


    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
}

