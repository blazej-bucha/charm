/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../src/prec.h"
#include "shc_touch_array_elements.h"
/* ------------------------------------------------------------------------- */






/* Check that "shcs->c" and "shcs->s" have access to the required number of
 * elements.  The bad thing is that if this condition is not satisfied,
 * a segfault error will appear on runtime, not leaving any room to report to
 * the user what kind of error was encountered.  Still, in this rare case,
 * a segfault error is still much better than no error at all. */
/* ------------------------------------------------------------------------- */
void shc_touch_array_elements(CHARM(shc) *shcs)
{
    if (shcs == NULL)
        return;


#if HAVE_MPI
    const size_t nchunk = shcs->local_nchunk;
    const unsigned long *orders = shcs->local_order;
#else
    const size_t nchunk = 1;
    const unsigned long orders[2] = {0, shcs->nmax};
#endif


    for (size_t k = 0; k < nchunk; k++)
    {
        for (unsigned long m = orders[2 * k]; m <= orders[2 * k + 1]; m++)
        {
            for (unsigned long n = m; n <= shcs->nmax; n++)
            {
                shcs->c[m][n - m] = PREC(1.0);
                shcs->s[m][n - m] = PREC(1.0);
            }
        }
    }


#if HAVE_MPI
    for (size_t j = 0; j < 2 * nchunk; j++)
        shcs->local_order[j] = 0;
#endif


    return;
}
/* ------------------------------------------------------------------------- */
