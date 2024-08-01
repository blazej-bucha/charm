/* This header file is not a part of API.
 *
 * Defines macros to compute sectorial, semi-sectorial and tesseral Legendre
 * functions.  We prefer macros over functions in this case to make sure the
 * code is always inlined in order to reduce the overhead.  The macros look
 * nasty, but significantly simplify things.
 *
 * */


#ifndef __LEG_FUNC_XNUM_H__
#define __LEG_FUNC_XNUM_H__


#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"


#ifdef __cplusplus
extern "C"
{
#endif


/* Macros to compute Legendre functions using "F"-numbers */
/* ------------------------------------------------------------------------- */
/* Tesseral Legendre functions */
#define PNM_RECURRENCE(x, y, pnm2, t, anms, bnms)                             \
        (pnm2) = SUB_R(MUL_R(MUL_R((anms), (t)), (x)),                        \
                       MUL_R((bnms), (y)));


/* First-order derivatives of tesseral Legendre functions */
#define DPNM_RECURRENCE(dpnm2, pnm1, pnm2, tu, u_rec, ns, enms)               \
        (dpnm2) = SUB_R(MUL_R(MUL_R((enms), (u_rec)), (pnm1)),                \
                        MUL_R(MUL_R((ns), (tu)), (pnm2)));


/* Second-order derivative of tesseral legendre functions */
#define DDPNM_RECURRENCE(ddpnm, dpnm, pnm, tu, u2_rec, m2s, nn1s)             \
        (ddpnm) = ADD_R(MUL_R((tu), (dpnm)),                                  \
                        MUL_R(SUB_R(MUL_R((m2s), (u2_rec)), (nn1s)), (pnm)));


/* Update the terms of a three-term recurrence */
#define RECURRENCE_NEXT_ITER(a0, a1, a2)                                      \
        {                                                                     \
            a0 = a1;                                                          \
            a1 = a2;                                                          \
        }
/* ------------------------------------------------------------------------- */






#ifdef SIMD
    /* Macros to compute sectorial, semisectorial and tesseral Legendre
     * functions using a SIMD code */
    /* --------------------------------------------------------------------- */
    /* Macro to compute sectorial Legendre functions. */
#   define PNM_SECTORIAL_XNUM_SIMD(x, ix, psmm1, ipsmm1, pnm0, BIG_r,         \
                                   zero_r, zero_ri, mone_ri, mask1, mask2,    \
                                   goto_label)                                \
                                                                              \
            (x) = LOAD_R(  &(psmm1));                                         \
           (ix) = LOAD_RI(&(ipsmm1));                                         \
                                                                              \
                                                                              \
           (mask1) = EQ_RI((ix), (zero_ri));                                  \
           (pnm0)  = BLEND_R((zero_r), (x), CAST_RI2R((mask1)));              \
           (mask2) = (mask1);                                                 \
           if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                               \
               goto goto_label;                                               \
                                                                              \
                                                                              \
           (mask1) = GT_RI((mone_ri), (ix));                                  \
           (pnm0)  = BLEND_R((pnm0), (zero_r), CAST_RI2R((mask1)));           \
           (mask2) = OR_MASK((mask1), (mask2));                               \
           if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                               \
               goto goto_label;                                               \
                                                                              \
                                                                              \
           (mask1) = EQ_RI((ix), (mone_ri));                                  \
           (pnm0)  = BLEND_R((pnm0), MUL_R((x), (BIGI_r)),                    \
                             CAST_RI2R((mask1)));                             \
           (mask2) = OR_MASK((mask1), (mask2));                               \
           if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                               \
               goto goto_label;                                               \
                                                                              \
                                                                              \
           (pnm0) = BLEND_R((pnm0), MUL_R((x), (BIG_r)),                      \
                            CAST_RI2R(ANDNOT_MASK((mask2))));                 \
                                                                              \
                                                                              \
    goto_label:






    /* Macro to compute semi-sectorial Legendre functions. */
#   define PNM_SEMISECTORIAL_XNUM_SIMD(x, y, ix, iy, w, t, anmmp1, pnm1,      \
                                       mask1, mask2, mask3,                   \
                                       zero_r, zero_ri, mone_ri,              \
                                       BIG_r, BIGS_r,  BIGI_r, goto_label)    \
                                                                              \
                                                                              \
            (y) =  (x);                                                       \
           (iy) = (ix);                                                       \
            (x) = MUL_R(MUL_R((anmmp1), (t)), (y));                           \
           (ix) = (iy);                                                       \
                                                                              \
                                                                              \
           (w)    = ABS_R((x));                                               \
           (mask3) = GE_R((w), (BIGS_r));                                     \
           if (MASK_TRUE_ANY(mask3))                                          \
           {                                                                  \
                (x) = BLEND_R((x), MUL_R((x), (BIGI_r)), (mask3));            \
               (ix) = BLEND_RI((ix), ADD_RI((ix), (one_ri)),                  \
                               CAST_R2RI((mask3)));                           \
           }                                                                  \
           (mask3) = LT_R((w), (BIGSI_r));                                    \
           if (MASK_TRUE_ANY(mask3))                                          \
           {                                                                  \
                (x) = BLEND_R((x), MUL_R((x), (BIG_r)), (mask3));             \
               (ix) = BLEND_RI((ix), SUB_RI((ix), (one_ri)),                  \
                               CAST_R2RI((mask3)));                           \
           }                                                                  \
                                                                              \
                                                                              \
           (mask1) = EQ_RI((ix), (zero_ri));                                  \
           (pnm1)  = BLEND_R((zero_r), (x), CAST_RI2R((mask1)));              \
           (mask2) = (mask1);                                                 \
           if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                               \
               goto goto_label;                                               \
                                                                              \
                                                                              \
           (mask1) = GT_RI((mone_ri), (ix));                                  \
           (pnm1)  = BLEND_R((pnm1), (zero_r), CAST_RI2R((mask1)));           \
           (mask2) = OR_MASK((mask1), (mask2));                               \
           if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                               \
               goto goto_label;                                               \
                                                                              \
                                                                              \
           (mask1) = EQ_RI((ix), (mone_ri));                                  \
           (pnm1)  = BLEND_R((pnm1), MUL_R((x), (BIGI_r)),                    \
                             CAST_RI2R((mask1)));                             \
           (mask2) = OR_MASK((mask1), (mask2));                               \
           if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                               \
               goto goto_label;                                               \
                                                                              \
                                                                              \
           (pnm1) = BLEND_R((pnm1), MUL_R((x), (BIG_r)),                      \
                            CAST_RI2R(ANDNOT_MASK((mask2))));                 \
                                                                              \
                                                                              \
    goto_label:






    /* Macro to compute tesseral Legendre functions. */
#   define PNM_TESSERAL_XNUM_SIMD(x, y, z, ix, iy, iz, ixy, w, t,             \
                                  anms, bnms, pnm2, tmp1_r, tmp2_r,           \
                                  mask1, mask2, mask3,                        \
                                  zero_r, zero_ri, one_ri,                    \
                                  BIG_r, BIGI_r, BIGS_r, BIGSI_r,             \
                                  goto_label1, goto_label2, ds)               \
                                                                              \
                                                                              \
           if ((ds))                                                          \
           {                                                                  \
               PNM_RECURRENCE(x, y, pnm2, t, anms, bnms);                     \
               RECURRENCE_NEXT_ITER(y, x, pnm2);                              \
           }                                                                  \
           else                                                               \
           {                                                                  \
               (ixy) = SUB_RI((ix), (iy));                                    \
                                                                              \
                                                                              \
               (mask1)  = EQ_RI((ixy), (zero_ri));                            \
               (tmp1_r) = BLEND_R((zero_r), (x), CAST_RI2R((mask1)));         \
               (tmp2_r) = BLEND_R((y), (y), CAST_RI2R((mask1)));              \
               (iz)     = BLEND_RI((zero_ri), (ix), (mask1));                 \
               (mask2)  = (mask1);                                            \
               if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                           \
                   goto goto_label1;                                          \
                                                                              \
                                                                              \
               (mask1)  = EQ_RI((ixy), (one_ri));                             \
               (tmp1_r) = BLEND_R((tmp1_r), (x), CAST_RI2R((mask1)));         \
               (tmp2_r) = BLEND_R((tmp2_r), MUL_R((y), (BIGI_r)),             \
                                  CAST_RI2R((mask1)));                        \
               (iz)     = BLEND_RI((iz), (ix), (mask1));                      \
               (mask2)  = OR_MASK((mask1), (mask2));                          \
               if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                           \
                   goto goto_label1;                                          \
                                                                              \
                                                                              \
               (mask1)  = EQ_RI((ixy), (mone_ri));                            \
               (tmp1_r) = BLEND_R((tmp1_r), MUL_R((x), (BIGI_r)),             \
                                  CAST_RI2R((mask1)));                        \
               (tmp2_r) = BLEND_R((tmp2_r), (y), CAST_RI2R((mask1)));         \
               (iz)     = BLEND_RI((iz), (iy), (mask1));                      \
               (mask2)  = OR_MASK((mask1), (mask2));                          \
               if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                           \
                   goto goto_label1;                                          \
                                                                              \
                                                                              \
               (mask1)  = GT_RI((ixy), (one_ri));                             \
               (tmp1_r) = BLEND_R((tmp1_r), (x), CAST_RI2R((mask1)));         \
               (tmp2_r) = BLEND_R((tmp2_r), (zero_r), CAST_RI2R((mask1)));    \
               (iz)     = BLEND_RI((iz), (ix), (mask1));                      \
               (mask2)  = OR_MASK((mask1), (mask2));                          \
               if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                           \
                   goto goto_label1;                                          \
                                                                              \
                                                                              \
               (iz) = BLEND_RI((iz), (iy), ANDNOT_MASK((mask2)));             \
                                                                              \
                                                                              \
    goto_label1:                                                              \
               (z) = SUB_R(MUL_R(MUL_R((anms), (t)), (tmp1_r)),               \
                           MUL_R((bnms), (tmp2_r)));                          \
                                                                              \
                                                                              \
               (w)    = ABS_R((z));                                           \
               (mask3) = GE_R((w), (BIGS_r));                                 \
               if (MASK_TRUE_ANY(mask3))                                      \
               {                                                              \
                    (z) = BLEND_R((z), MUL_R((z), (BIGI_r)), (mask3));        \
                   (iz) = BLEND_RI((iz), ADD_RI((iz), (one_ri)),              \
                                   CAST_R2RI((mask3)));                       \
               }                                                              \
               (mask3) = LT_R((w), (BIGSI_r));                                \
               if (MASK_TRUE_ANY(mask3))                                      \
               {                                                              \
                    (z) = BLEND_R((z), MUL_R((z), (BIG_r)), (mask3));         \
                   (iz) = BLEND_RI((iz), SUB_RI((iz), (one_ri)),              \
                                   CAST_R2RI((mask3)));                       \
               }                                                              \
                                                                              \
                                                                              \
                (y) =  (x);                                                   \
               (iy) = (ix);                                                   \
                (x) =  (z);                                                   \
               (ix) = (iz);                                                   \
                                                                              \
                                                                              \
               (mask1) = EQ_RI((iz), (zero_ri));                              \
               (pnm2)  = BLEND_R((zero_r), (z), CAST_RI2R((mask1)));          \
               (mask2) = (mask1);                                             \
               if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                           \
                   goto goto_label2;                                          \
                                                                              \
                                                                              \
               (mask1) = GT_RI((mone_ri), (iz));                              \
               (pnm2)  = BLEND_R((pnm2), (zero_r), CAST_RI2R((mask1)));       \
               (mask2) = OR_MASK((mask1), (mask2));                           \
               if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                           \
                   goto goto_label2;                                          \
                                                                              \
                                                                              \
               (mask1) = EQ_RI((iz), (mone_ri));                              \
               (pnm2)  = BLEND_R((pnm2), MUL_R((z), (BIGI_r)),                \
                                 CAST_RI2R((mask1)));                         \
               (mask2) = OR_MASK((mask1), (mask2));                           \
               if (MASK_TRUE_ALL(CAST_RI2R(mask2)))                           \
                   goto goto_label2;                                          \
                                                                              \
                                                                              \
               (pnm2) = BLEND_R((pnm2), MUL_R((z), (BIG_r)),                  \
                                CAST_RI2R(ANDNOT_MASK((mask2))));             \
                                                                              \
                                                                              \
    goto_label2:                                                              \
               if (MASK_TRUE_ALL(CAST_RI2R(EQ_RI((ix), (zero_ri)))) &&        \
                   MASK_TRUE_ALL(CAST_RI2R(EQ_RI((iy), (zero_ri)))))          \
                   /* Excellent, dynamical switching can be applied from */   \
                   /* now on for all degrees, so no X-numbers! */             \
                   /* IMPORTANTLY, this may not work with other */            \
                   /* normalization schemes that may be added in the */       \
                   /* future. */                                              \
                   (ds) = 1;                                                  \
        }
#endif
/* ------------------------------------------------------------------------- */






/* Macros to compute sectorial, semisectorial and tesseral Legendre functions
 * using a scalar code */
/* ------------------------------------------------------------------------- */
/* Macro to compute sectorial Legendre functions. */
#define PNM_SECTORIAL_XNUM(x, ix, psmm1, ipsmm1, pnm0)                        \
                                                                              \
                                                                              \
         (x) = (psmm1);                                                       \
        (ix) = (ipsmm1);                                                      \
                                                                              \
                                                                              \
        if ((ix) == 0)                                                        \
            (pnm0) = (x);                                                     \
        else if ((ix) < -1)                                                   \
            (pnm0) = PREC(0.0);                                               \
        else if ((ix) < 0)                                                    \
            (pnm0) = (x) * BIGI;                                              \
        else                                                                  \
            (pnm0) = (x) * BIG;






/* Macro to compute semi-sectorial Legendre functions. */
#define PNM_SEMISECTORIAL_XNUM(x, y, ix, iy, w, t, anm, pnm1)                 \
                                                                              \
                                                                              \
         (y) =  (x);                                                          \
        (iy) = (ix);                                                          \
         (x) = ((anm) * (t)) * (y);                                           \
        (ix) = (iy);                                                          \
                                                                              \
                                                                              \
        (w) = FABS((x));                                                      \
        if ((w) >= BIGS)                                                      \
        {                                                                     \
             (x) *= BIGI;                                                     \
            (ix) += 1;                                                        \
        }                                                                     \
        else if ((w) < BIGSI)                                                 \
        {                                                                     \
             (x) *= BIG;                                                      \
            (ix) -= 1;                                                        \
        }                                                                     \
                                                                              \
                                                                              \
        if ((ix) == 0)                                                        \
            (pnm1) = (x);                                                     \
        else if ((ix) < -1)                                                   \
            (pnm1) = PREC(0.0);                                               \
        else if ((ix) < 0)                                                    \
            (pnm1) = (x) * BIGI;                                              \
        else                                                                  \
            (pnm1) = (x) * BIG;






/* Macro to compute tesseral Legendre functions.  Typical usage is inside the
 * innermost loop over harmonic degrees.
 *
 *      PNM_TESSERAL_XNUM(x, y, z, ix, iy, iz, ixy, w, t, anm[n], bnm[n], pnm2,
 *                        continue)
 *
 * The computed tesseral Legendre function will be stored in "pnm2".
 *
 * */
#define PNM_TESSERAL_XNUM(x, y, z, ix, iy, iz, ixy, w, t,                     \
                          anm, bnm, pnm2, ds)                                 \
                                                                              \
                                                                              \
        if ((ds))                                                             \
        {                                                                     \
            PNM_RECURRENCE(x, y, pnm2, t, anm, bnm);                          \
            RECURRENCE_NEXT_ITER(y, x, pnm2);                                 \
        }                                                                     \
        else                                                                  \
        {                                                                     \
            (ixy) = (ix) - (iy);                                              \
                                                                              \
                                                                              \
            if ((ixy) == 0)                                                   \
            {                                                                 \
                 (z) = ((anm) * (t)) * (x) - (bnm) * (y);                     \
                (iz) = (ix);                                                  \
            }                                                                 \
            else if ((ixy) == 1)                                              \
            {                                                                 \
                 (z) = ((anm) * (t)) * (x) - (bnm) * ((y) * BIGI);            \
                (iz) = (ix);                                                  \
            }                                                                 \
            else if ((ixy) == -1)                                             \
            {                                                                 \
                 (z) = ((anm) * (t)) * ((x) * BIGI) - (bnm) * (y);            \
                (iz) = (iy);                                                  \
            }                                                                 \
            else if ((ixy) > 1)                                               \
            {                                                                 \
                 (z) = ((anm) * (t)) * (x);                                   \
                (iz) = (ix);                                                  \
            }                                                                 \
            else                                                              \
            {                                                                 \
                 (z) = -(bnm) * (y);                                          \
                (iz) = (iy);                                                  \
            }                                                                 \
                                                                              \
                                                                              \
            (w) = FABS((z));                                                  \
            if ((w) >= BIGS)                                                  \
            {                                                                 \
                 (z) *= BIGI;                                                 \
                (iz) += 1;                                                    \
            }                                                                 \
            else if ((w) < BIGSI)                                             \
            {                                                                 \
                 (z) *= BIG;                                                  \
                (iz) -= 1;                                                    \
            }                                                                 \
                                                                              \
                                                                              \
             (y) =  (x);                                                      \
            (iy) = (ix);                                                      \
             (x) =  (z);                                                      \
            (ix) = (iz);                                                      \
                                                                              \
                                                                              \
            if ((iz) == 0)                                                    \
                (pnm2) = (z);                                                 \
            else if ((iz) < -1)                                               \
                (pnm2) = PREC(0.0);                                           \
            else if ((iz) < 0)                                                \
                (pnm2) = (z) * BIGI;                                          \
            else                                                              \
                (pnm2) = (z) * BIG;                                           \
                                                                              \
                                                                              \
            if (((ix) == 0) && ((iy) == 0))                                   \
                /* Excellent, dynamical switching can be applied from */      \
                /* now on for all degrees, so no X-numbers! */                \
                /* IMPORTANTLY, this may not work with other */               \
                /* normalization schemes that may be added in the */          \
                /* future. */                                                 \
                (ds) = 1;                                                     \
        }
/* ------------------------------------------------------------------------- */


#ifdef __cplusplus
}
#endif


#endif

