/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






size_t CHARM(shs_cell_isurf_offset)(unsigned long d1, unsigned long d2,
                                    unsigned long d3, unsigned long d4,
                                    size_t nd2, size_t nd3, size_t nd4)
/*
 * ============================================================================
 *
 * DESCRIPTION: Function to index elements of a 4D array, "d1", "d2", "d3",
 *              "d4", that has dimmensions "nd1", "nd2", "nd3", "nd4".
 *
 * ============================================================================
 *
 * */
{
    return ((d1 * nd2 + d2) * nd3 + d3) * nd4 + d4;
}
