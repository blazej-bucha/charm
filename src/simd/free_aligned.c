/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include "simd.h"
#include "free_aligned.h"
/* ------------------------------------------------------------------------- */






/* Frees the memory dynamically allocated by "CHARM(malloc_aligned)" or
 * "CHARM(calloc_aligned)". */
void CHARM(free_aligned)(void *ptr)
{
#ifndef SIMD
    free(ptr);
    return;
#elif HAVE_POSIX_MEMALIGN
    free(ptr);
    return;
#elif HAVE_ALIGNED_ALLOC
    free(ptr);
    return;
#elif HAVE_MM_MALLOC_H
    _mm_free(ptr);
    return;
#else
#   error "Couldn't find any system function to free an aligned block of memory."
#endif
}
