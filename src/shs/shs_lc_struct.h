/* This header file is not a part of API. */


#ifndef __SHS_LC_STRUCT_H__
#define __SHS_LC_STRUCT_H__


#include "../prec.h"
#include "../simd/simd.h"


#undef LC_BLOCKS
#define LC_BLOCKS (24)


/* Structure to store the lumped coefficients. */
typedef struct
{
#if HAVE_MPI
    /* In this case, the memory will be allocated and freed dynamically */
    REAL_SIMD *_all;
#else
    /* In this case, "_all" is allocated statically.
     *
     * Importantly, with SIMD instructions it is absolutely crucial for the
     * "_all" member to be the first member of the structure, because it is of
     * "REAL_SIMD" data type, so its memory must be aligned.  We achieve this
     * by considering that the C standard guarantees that the pointers to
     * a structure and to its first member are equal.  So all we need to do is
     * to allocate aligned memory for the whole structure, see "shs_lc_init.c"
     * for details. */
    REAL_SIMD _all[LC_BLOCKS * SIMD_BLOCK_S];
#endif

    REAL_SIMD *a;     /* "pnm" */
    REAL_SIMD *b;     /* "pnm" */
    REAL_SIMD *a2;    /* "pnm" */
    REAL_SIMD *b2;    /* "pnm" */


    REAL_SIMD *ar;    /* "(n + 1) * pnm" */
    REAL_SIMD *ap;    /* "dpnm" */
    REAL_SIMD *arr;   /* "(n + 1) * (n + 2) * pnm" */
    REAL_SIMD *arp;   /* "(n + 1) * dpnm" */
    REAL_SIMD *app;   /* "ddpnm" */


    REAL_SIMD *br;    /* "(n + 1) * pnm" */
    REAL_SIMD *bp;    /* "dpnm" */
    REAL_SIMD *brr;   /* "(n + 1) * (n + 2) * pnm" */
    REAL_SIMD *brp;   /* "(n + 1) * dpnm" */
    REAL_SIMD *bpp;   /* "ddpnm" */


    REAL_SIMD *ar2;   /* "(n + 1) * pnm" */
    REAL_SIMD *ap2;   /* "dpnm" */
    REAL_SIMD *arr2;  /* "(n + 1) * (n + 2) * pnm" */
    REAL_SIMD *arp2;  /* "(n + 1) * dpnm" */
    REAL_SIMD *app2;  /* "ddpnm" */


    REAL_SIMD *br2;   /* "(n + 1) * pnm" */
    REAL_SIMD *bp2;   /* "dpnm" */
    REAL_SIMD *brr2;  /* "(n + 1) * (n + 2) * pnm" */
    REAL_SIMD *brp2;  /* "(n + 1) * dpnm" */
    REAL_SIMD *bpp2;  /* "ddpnm" */


    /* The shs point kernels do not use the error structure, so this variable
     * is used to indicate an error in the kernels back to the caller. */
    _Bool error;
} CHARM(lc);


#endif
