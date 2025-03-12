/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Returns index of the "d1", "d2", "d3", "d4" element of a 4D array of
 * dimmensions "nd1", "nd2", "nd3", "nd4".  */
size_t CHARM(misc_idx_4d)(unsigned long d1,
                          unsigned long d2,
                          unsigned long d3,
                          unsigned long d4,
                          size_t nd2,
                          size_t nd3,
                          size_t nd4)
{
    return ((d1 * nd2 + d2) * nd3 + d3) * nd4 + d4;
}
