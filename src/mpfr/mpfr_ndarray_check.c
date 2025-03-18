/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include <stdio.h>
#include <stdarg.h>
#include <mpfr.h>
#include "mpfr_ndarray.h"
/* ------------------------------------------------------------------------- */






/* Internal function to check whether an "mpfr_ndarray" "x" is of at least
 * "ndim" dimensions and whether the number of elements stored in these
 * dimensions is sufficiently large.
 *
 * IMPORTANT NOTE: See "mpfr_ndarray_malloc" on explicit casting of the
 * variadic input parameters. */
int CHARM(mpfr_ndarray_check)(const mpfr_ndarray *x,
                              size_t ndim,
                              ...)
{
    if (x == NULL)
        return 1;


    /* If "ndim" is set to zero, this means do not check the dimensions of the
     * ndarray */
    if (ndim == 0)
        return 0;


    va_list ap;
    va_start(ap, ndim);
    size_t n;
    int ret = 0;
    for (size_t d = 0; d < ndim; d++)
    {
        n = va_arg(ap, size_t);


        /* Zero "n" means do not check that particular dimension of the ndarray
         * */
        if ((n > 0) && (x->shape[d] < n))
        {
            ret = 2;
            goto EXIT;
        }
    }


EXIT:
    va_end(ap);
    return ret;
}
