/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#if HAVE__ALIGNED_MALLOC
#   include <malloc.h>
#endif
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
#elif HAVE_POSIX_MEMALIGN
    free(ptr);
#elif HAVE_ALIGNED_ALLOC
    free(ptr);
#elif HAVE_MM_MALLOC_H
    _mm_free(ptr);
#elif HAVE__ALIGNED_MALLOC
    _aligned_free(ptr);
#else
#   error "Couldn't find any system function to free an aligned block of memory."
#endif
    return;
}
