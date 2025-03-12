/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <mpfr.h>
#include "mpfr_flush_unreleased_memory.h"
#include "mpfr_ndarray.h"
#include "mpfr_check_bits.h"
/* ------------------------------------------------------------------------- */







/* Internal variadic function to allocate a "mpfr_ndarray" of "ndim" dimensions
 * using "NBITS" for "mpfr_t".  After "ndim" follow numbers of elements to be
 * stored in each dimension.  For instance, to allocate an array of two
 * dimensions with "5" and "6" elements, use
 *
 *      CHARM(mpfr_ndarray_malloc)(NBITS, 2, 5, 6);
 *
 * */
mpfr_ndarray *CHARM(mpfr_ndarray_malloc)(mpfr_prec_t NBITS,
                                         size_t ndim,
                                         ...)
{
    /* Check inputs */
    /* --------------------------------------------------------------------- */
    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
        return NULL;


    CHARM(mpfr_check_bits)(NBITS, err);
    _Bool NBITS_correct = CHARM(err_isempty)(err) ? 1 : 0;
    CHARM(err_free)(err);
    if (!NBITS_correct)
        return NULL;


    if (ndim < 1)
        return NULL;
    /* --------------------------------------------------------------------- */


    /* Allocate memory for the "mpfr_ndarray" data type */
    /* --------------------------------------------------------------------- */
    mpfr_ndarray *x = malloc(sizeof(mpfr_ndarray));
    if (x == NULL)
        return NULL;


    x->data  = NULL;
    x->shape = NULL;
    /* --------------------------------------------------------------------- */


    /* Set the "ndim" member */
    /* --------------------------------------------------------------------- */
    x->ndim = ndim;
    /* --------------------------------------------------------------------- */


    /* Set the "shape" member */
    /* --------------------------------------------------------------------- */
    va_list ap;
    va_start(ap, ndim);


    x->shape = (size_t *)malloc(x->ndim * sizeof(size_t));
    if (x->shape == NULL)
        goto FAILURE;


    for (size_t d = 0; d < x->ndim; d++)
    {
        x->shape[d] = va_arg(ap, size_t);
        if (x->shape[d] < 1)
            goto FAILURE;
    }
    /* --------------------------------------------------------------------- */


    /* Allocate data memory and initialize all array elements */
    /* --------------------------------------------------------------------- */
    /* Get the total number of elements that need to be allocated */
    x->size = 1;
    for (size_t d = 0; d < x->ndim; d++)
        x->size *= x->shape[d];


    /* Allocate the memory for all array elements to be stored in the "data"
     * member */
    x->data = (mpfr_t *)malloc(x->size * sizeof(mpfr_t));
    if (x->data == NULL)
        goto FAILURE;


    /* Initialize all array elements */
    for (size_t i = 0; i < x->size; i++)
        mpfr_init2(x->data[i], NBITS);
    /* --------------------------------------------------------------------- */


    /* Set the owner member to "1" to signalize that "mpfr_ndarray_free" is
     * allowed to deallocate the memory associated with "x->data" if asked to
     * do so. */
    /* --------------------------------------------------------------------- */
    x->owner = 1;
    /* --------------------------------------------------------------------- */


    /* Normal exit */
    /* --------------------------------------------------------------------- */
EXIT:
    va_end(ap);
    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return x;
    /* --------------------------------------------------------------------- */


    /* Premature exit */
    /* --------------------------------------------------------------------- */
FAILURE:
    free(x->data);
    free(x->shape);
    free(x);
    x = NULL;


    goto EXIT;
    /* --------------------------------------------------------------------- */
}

