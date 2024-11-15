/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "shc_local_ncs.h"
/* ------------------------------------------------------------------------- */






/* Returns the total number of coefficients to be stored locally based on the
 * chunk specification.  Must be called only if "nmax", "nchunk" and "order"
 * were already checked by "shc_check_chunk_orders". */
size_t CHARM(shc_local_ncs)(unsigned long nmax,
                            size_t nchunk,
                            const unsigned long *order)
{
    unsigned long nmax1 = (size_t)(nmax) + 1;
    size_t ncs = 0;
    for (size_t j = 0; j < nchunk; j++)
        for (unsigned long m = order[2 * j]; m <= order[2 * j + 1]; m++)
            ncs += nmax1 - m;


    return ncs;
}
