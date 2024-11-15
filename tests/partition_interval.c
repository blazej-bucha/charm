/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdint.h>
#include <limits.h>
#include "../src/prec.h"
#include "partition_interval.h"
/* ------------------------------------------------------------------------- */






/* Computes respectively the minimum and maximum values "out_start" and
 * "out_stop" of the "i"th partition, "i = 0, 1, ..., count - 1" of the
 * "[start, stop]" interval.  */
/* ------------------------------------------------------------------------- */
#undef TEMPLATE
#define TEMPLATE(FUNC, TYPE, MAX_VAL) void FUNC(TYPE start,                   \
                                                TYPE stop,                    \
                                                TYPE count,                   \
                                                TYPE i,                       \
                                                TYPE *out_start,              \
                                                TYPE *out_stop)               \
{                                                                             \
    if (i >= count)                                                           \
    {                                                                         \
        *out_start = *out_stop = MAX_VAL;                                     \
        return;                                                               \
    }                                                                         \
                                                                              \
                                                                              \
    unsigned long chunk = (stop - start) / count;                             \
    *out_start = (i == 0) ? start : i * chunk + start + 1;                    \
    *out_stop = (i == (count - 1)) ? stop : (i + 1) * chunk + start;          \
                                                                              \
                                                                              \
    if (*out_start > *out_stop)                                               \
    {                                                                         \
        *out_start = *out_stop = MAX_VAL;                                     \
        return;                                                               \
    }                                                                         \
                                                                              \
                                                                              \
    return;                                                                   \
}


TEMPLATE(partition_interval_ulong, unsigned long, ULONG_MAX)
TEMPLATE(partition_interval_size_t, size_t, SIZE_MAX)


#undef TEMPLATE
/* ------------------------------------------------------------------------- */

