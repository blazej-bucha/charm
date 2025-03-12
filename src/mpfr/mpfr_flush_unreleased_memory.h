/* This header file is not a part of API. */


#ifndef __MPFR_FLUSH_UNRELEASED_MEMORY_H__
#define __MPFR_FLUSH_UNRELEASED_MEMORY_H__


#undef FLUSH_UNRELEASED_MEMORY
#if defined(__GLIBC__) && HAVE_MALLOC_H && HAVE_MALLOC_TRIM
    /* If we got here, it means we linked to "glibc", we have "malloc.h" and we
     * have the "malloc_trim" function. */
#   include <malloc.h>


    /* "glibc" with the default settings for "malloc" allocates memory for
     * arrays of "mpfr_t" data type in such a way that not all of the allocated
     * memory is properly freed (returned back to the system) when calling
     * "free".  In fact, if the number of "mpfr_t" array elements is large,
     * this may cause GBs of memory to be not freed properly.  This happens
     * despite the code to allocate and deallocate the arrays does not have any
     * memory leaks.  This was confirmed by "valgrind" but also by linking
     * against the musl's "libc", which works much better out-of-the-box.
     * A workaround is to call "malloc_trim" after all the commands to release
     * the memory associated with the "mpfr_t" arrays have been executed.
     *
     * The issue happens with arrays of the "mpfr_t" data type, but MPFR is not
     * to be blamed for this.  This is indeed caused by "glibc" and the problem
     * seems to be in that calling "mpfr_init2" in a loop means that MPFR
     * dynamically allocates small amounts of memory (or large, depending on
     * "mpfr_prec_t") with each iteration.  These many dynamical allocations of
     * small memory chunks then cause the troubles when using the default
     * "malloc" parameters of "glibc" (see "mallopt").
     *
     * Importantly, "malloc_trim" is a GNU extension, so we must check whether
     * "__GLIBC__" is defined.  If "__GLIBC__" is not defined, this is OK,
     * given that it is the glibc library, which causes the issue in the first
     * place.
     *
     * In practice, the "FLUSH_UNRELEASED_MEMORY" macro should be called
     * whenever deallocating arrays of "mpfr_t". */
#   define FLUSH_UNRELEASED_MEMORY malloc_trim(0);
#else
#   define FLUSH_UNRELEASED_MEMORY
#endif


#endif

