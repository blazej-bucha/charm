/* This header file is not a part of API.
 *
 * The macros allocate arrays either statically or dynamically, depending on
 * whether or not we are compiling with MPI support.
 *
 * */


#ifndef __MISC_SD_CALLOC_H__
#define __MISC_SD_CALLOC_H__


#include <config.h>
#include "../prec.h"
#include "../src/simd/simd.h"
#include "../src/simd/calloc_aligned.h"
#include "../src/simd/free_aligned.h"


/* Initializations */
/* ------------------------------------------------------------------------- */
#undef MISC_SD_CALLOC_REAL_SIMD_INIT
#if HAVE_MPI
#   define MISC_SD_CALLOC_REAL_SIMD_INIT(x) REAL_SIMD *(x) = NULL;
#else
#   define MISC_SD_CALLOC_REAL_SIMD_INIT(x)
#endif


#undef MISC_SD_CALLOC_RI_SIMD_INIT
#if HAVE_MPI
#   define MISC_SD_CALLOC_RI_SIMD_INIT(x) RI_SIMD *(x) = NULL;
#else
#   define MISC_SD_CALLOC_RI_SIMD_INIT(x)
#endif


#undef MISC_SD_CALLOC__BOOL_INIT
#if HAVE_MPI
#   define MISC_SD_CALLOC__BOOL_INIT(x) _Bool *(x) = NULL;
#else
#   define MISC_SD_CALLOC__BOOL_INIT(x)
#endif
/* ------------------------------------------------------------------------- */


/* Allocation using the error structure */
/* ------------------------------------------------------------------------- */
#undef MISC_SD_CALLOC_REAL_SIMD_ERR
#if HAVE_MPI
#   define MISC_SD_CALLOC_REAL_SIMD_ERR(x, size_dyn, size_stat, err,          \
                                        goto_label)                           \
    (x) = (REAL_SIMD *)CHARM(calloc_aligned)(SIMD_MEMALIGN, (size_dyn),       \
                                             sizeof(REAL_SIMD));              \
    if ((x) == NULL)                                                          \
    {                                                                         \
        CHARM(err_set)((err), __FILE__, __LINE__, __func__, CHARM_EMEM,       \
                       CHARM_ERR_MALLOC_FAILURE);                             \
        goto goto_label;                                                      \
    }
#else
#   define MISC_SD_CALLOC_REAL_SIMD_ERR(x, size_dyn, size_stat, err,          \
                                        goto_label)                           \
        REAL_SIMD (x)[(size_stat)];
#endif


#undef MISC_SD_CALLOC_RI_SIMD_ERR
#if HAVE_MPI
#   define MISC_SD_CALLOC_RI_SIMD_ERR(x, size_dyn, size_stat, err,            \
                                      goto_label)                             \
    (x) = (RI_SIMD *)CHARM(calloc_aligned)(SIMD_MEMALIGN, (size_dyn),         \
                                           sizeof(RI_SIMD));                  \
    if ((x) == NULL)                                                          \
    {                                                                         \
        CHARM(err_set)((err), __FILE__, __LINE__, __func__, CHARM_EMEM,       \
                       CHARM_ERR_MALLOC_FAILURE);                             \
        goto goto_label;                                                      \
    }
#else
#   define MISC_SD_CALLOC_RI_SIMD_ERR(x, size_dyn, size_stat, err, goto_label)\
        RI_SIMD (x)[(size_stat)];
#endif


#undef MISC_SD_CALLOC__BOOL_ERR
#if HAVE_MPI
#   define MISC_SD_CALLOC__BOOL_ERR(x, size_dyn, size_stat, err, goto_label)  \
    (x) = (_Bool *)CHARM(calloc_aligned)(SIMD_MEMALIGN, (size_dyn),           \
                                         sizeof(RI_SIMD));                    \
    if ((x) == NULL)                                                          \
    {                                                                         \
        CHARM(err_set)((err), __FILE__, __LINE__, __func__, CHARM_EMEM,       \
                       CHARM_ERR_MALLOC_FAILURE);                             \
        goto goto_label;                                                      \
    }
#else
#   define MISC_SD_CALLOC__BOOL_ERR(x, size_dyn, size_stat, err, goto_label)  \
        _Bool (x)[(size_stat)];
#endif
/* ------------------------------------------------------------------------- */


/* Allocation using the error structure */
/* ------------------------------------------------------------------------- */
#undef MISC_SD_CALLOC_REAL_SIMD_E
#if HAVE_MPI
#   define MISC_SD_CALLOC_REAL_SIMD_E(x, size_dyn, size_stat, e,              \
                                          goto_label)                         \
    (x) = (REAL_SIMD *)CHARM(calloc_aligned)(SIMD_MEMALIGN, (size_dyn),       \
                                             sizeof(REAL_SIMD));              \
    if ((x) == NULL)                                                          \
    {                                                                         \
        (e) = 1;                                                              \
        goto goto_label;                                                      \
    }
#else
#   define MISC_SD_CALLOC_REAL_SIMD_E MISC_SD_CALLOC_REAL_SIMD_ERR
#endif


#undef MISC_SD_CALLOC_RI_SIMD_E
#if HAVE_MPI
#   define MISC_SD_CALLOC_RI_SIMD_E(x, size_dyn, size_stat, e, goto_label)    \
    (x) = (RI_SIMD *)CHARM(calloc_aligned)(SIMD_MEMALIGN, (size_dyn),         \
                                           sizeof(RI_SIMD));                  \
    if ((x) == NULL)                                                          \
    {                                                                         \
        (e) = 1;                                                              \
        goto goto_label;                                                      \
    }
#else
#   define MISC_SD_CALLOC_RI_SIMD_E MISC_SD_CALLOC_RI_SIMD_ERR
#endif


#undef MISC_SD_CALLOC__BOOL_E
#if HAVE_MPI
#   define MISC_SD_CALLOC__BOOL_E(x, size_dyn, size_stat, e, goto_label)      \
    (x) = (_Bool *)CHARM(calloc_aligned)(SIMD_MEMALIGN, (size_dyn),           \
                                         sizeof(RI_SIMD));                    \
    if ((x) == NULL)                                                          \
    {                                                                         \
        (e) = 1;                                                              \
        goto goto_label;                                                      \
    }
#else
#   define MISC_SD_CALLOC__BOOL_E MISC_SD_CALLOC__BOOL_ERR
#endif
/* ------------------------------------------------------------------------- */


/* Freeing the memory.  The same macro is used for all data types */
/* ------------------------------------------------------------------------- */
#undef MISC_SD_FREE
#if HAVE_MPI
#   define MISC_SD_FREE(x) CHARM(free_aligned)(x);
#else
#   define MISC_SD_FREE(x)
#endif
/* ------------------------------------------------------------------------- */


#endif
