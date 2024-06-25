/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#if HAVE__ALIGNED_MALLOC
#   include <malloc.h>
#endif
#include "../prec.h"
#include "simd.h"
#include "malloc_aligned.h"
/* ------------------------------------------------------------------------- */






/* Dynamically allocates an aligned block of memory.  The memory must be freed
 * with "CHARM(free_aligned)".
 *
 * * If compiling without the SIMD support, used is the standard "malloc"
 *   function with its alignment.  This is to make sure that CHarm compiles
 *   even if the special memory allocation functions (such as "posix_memalign")
 *   are not available on the host system (given that they are not needed in
 *   this case).
 *
 * * If compiling with the SIMD support, used will be the first system function
 *   to allocate an aligned block of memory from the following list that is
 *   found on the host: "posix_memalign", "aligned_alloc", "_mm_malloc",
 *   "_aligned_malloc".
 *
 * */
void *CHARM(malloc_aligned)(size_t alignment, size_t size)
{
#ifndef SIMD
    /* If compiling without the SIMD support, the function simply calls
     * "malloc".  The "alignment" variable must be set to "0" which merely
     * means the alignment will be the same as that of "malloc".  If
     * "alignment" is not zero, "NULL" is returned to make sure the caller does
     * not attempt to specify any alignment. */
    if (alignment != 0)
        return NULL;
    return malloc(size);
#elif HAVE_POSIX_MEMALIGN
    void *p;
    if (posix_memalign(&p, alignment, size))
        return NULL;
    return p;
#elif HAVE_ALIGNED_ALLOC
    /* With "aligned_alloc", the number of bytes to be allocated must be an
     * integer multiple of "alignment".  Therefore, we get at first the
     * smallest integer that is a multiple of "alignment" and is not smaller
     * than "size".  Then, an alligned memory is allocated with
     * "aligned_alloc". */
    size_t size_multiple = (size_t)ceil((double)size / (double)alignment) *
                                        alignment;
    return aligned_alloc(alignment, size_multiple);
#elif HAVE_MM_MALLOC_H
    return _mm_malloc(size, alignment);
#elif HAVE__ALIGNED_MALLOC
    return _aligned_malloc(size, alignment);
#else
#   error "Couldn't find any system function to allocate an aligned block of memory."
#endif
}
