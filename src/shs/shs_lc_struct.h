/* This header file is not a part of API. */


#ifndef __SHS_LC_STRUCT_H__
#define __SHS_LC_STRUCT_H__


#include "../prec.h"
#include "../simd/simd.h"


#undef LC_BLOCKS
#define LC_BLOCKS (24)


/* Structure to store the lumped coefficients.  Do not change the order of the
 * arrays, as there is some logic behind it to offer reasonable caching. */
typedef struct
{
    REAL_SIMD _all[LC_BLOCKS * SIMD_BLOCK_S];

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
} CHARM(lc);


#endif
