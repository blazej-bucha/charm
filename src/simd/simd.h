/*
 * This file is not a part of the API.
 *
 * Defines macros, functions, symbolic constants, etc., to work with SIMD.
 *
 * */


#ifndef __SIMD_H__
#define __SIMD_H__


/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <string.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* At first, let's check if one type of SIMD instruction only is defined in
 * "config.h".  The configure script does not allow this, but if CHarm is
 * compiled without the autotools, this might perhaps happen. */
#if (HAVE_AVX  && HAVE_AVX2) || \
    (HAVE_AVX2 && HAVE_AVX512F) || \
    (HAVE_AVX  && HAVE_AVX512F)
#   error "One type of SIMD instructions only can be defined in config.h."
#endif






#undef SIMD
#undef SIMD_SIZE
#undef SIMD_BLOCK
#undef SIMD_MEMALIGN
#undef SIMD_TRUE
#undef REAL_SIMD
#undef INT_SIMD
#undef RI_SIMD
#undef MASK_SIMD
#undef MASK2_SIMD
#undef P
#undef PF
#undef PINT
#undef PINT2
#undef PF_MASK
#undef MUL_R
#undef DIV_R
#undef ADD_R
#undef ADD_RI
#undef SUB_R
#undef SUB_RI
#undef NEG_R_INIT
#undef NEG_R
#undef SIGNBIT
#undef SET1_R
#undef SET1_RI
#undef LOAD_R
#undef LOADU_R
#undef SUM_R
#undef ABS_R
#undef ABS_R_INIT
#undef NONSIGNBITS
#undef SIMD_MULTIPLE
#undef CAST_RI2R
#undef CAST_R2RI
#undef SET_ZERO_R
#undef SET_ZERO_I
#undef SET_ZERO_RI
#undef STORE_R
#undef STOREU_R
#undef EQ_R
#undef GE_R
#undef LT_R
#undef EQ_RI
#undef GT_RI
#undef OR_R
#undef AND_R
#undef OR_MASK
#undef ANDNOT_MASK
#undef LOAD_RI
#undef BLEND_R
#undef BLEND_RI
#undef MOVEMASK






#if HAVE_AVX || HAVE_AVX2 || HAVE_AVX512F
    /* If "SIMD" is defined, CHarm is being compiled with SIMD support. */
#   define SIMD
#endif






#ifdef SIMD


#   if CHARM_QUAD
#       error "Cannot compile with SIMD instructions in quadruple precision."
#   endif


    /* Header files */
    /* --------------------------------------------------------------------- */
#   include <immintrin.h>
#   include <math.h>
    /* --------------------------------------------------------------------- */





#   if HAVE_AVX || HAVE_AVX2

#       if CHARM_FLOAT /* Single precision */

#           define SIMD_SIZE     8
#           define SIMD_TRUE     0xFF
#           define PF(x)         _mm256_ ## x ## _ps
#           define PINT(x)       _mm256_ ## x ## _epi32
#           define PINT2         PINT
#           define REAL_SIMD     __m256

#       else /* Double precision */

#           define SIMD_SIZE     4
#           define SIMD_TRUE     0xF
#           define PF(x)         _mm256_ ## x ## _pd
#           define PINT(x)       _mm256_ ## x ## _epi64
#           define PINT2(x)      _mm256_ ## x ## _epi64x
#           define REAL_SIMD     __m256d

#       endif


#       define SIMD_MEMALIGN     32
#       define P(x)              _mm256_ ## x
#       define INT_SIMD          __m256i


#       if HAVE_AVX
#           define RI_SIMD       REAL_SIMD
#       elif HAVE_AVX2
#           define RI_SIMD       INT_SIMD
#       endif


#       define MASK_SIMD         RI_SIMD
#       define MASK2_SIMD        REAL_SIMD

#   elif HAVE_AVX512F

#       if CHARM_FLOAT /* Single precision */

#           define SIMD_SIZE     16
#           define SIMD_TRUE     0xFFFF
#           define PF(x)         _mm512_ ## x ## _ps
#           define PF_MASK(x)    _mm512_ ## x ## _ps_mask
#           define PINT(x)       _mm512_ ## x ## _epi32
#           define PI_MASK(x)    _mm512_ ## x ## _epi32_mask
#           define REAL_SIMD     __m512
#           define MASK_SIMD     __mmask16
#           define MASK2_SIMD    MASK_SIMD

#       else /* Double precision */

#           define SIMD_SIZE     8
#           define SIMD_TRUE     0xFF
#           define PF(x)         _mm512_ ## x ## _pd
#           define PF_MASK(x)    _mm512_ ## x ## _pd_mask
#           define PINT(x)       _mm512_ ## x ## _epi64
#           define PI_MASK(x)    _mm512_ ## x ## _epi64_mask
#           define REAL_SIMD     __m512d
#           define MASK_SIMD     __mmask8
#           define MASK2_SIMD    MASK_SIMD

#       endif


#       define SIMD_MEMALIGN     64
#       define PINT2             PINT
#       define P(x)              _mm512_ ## x
#       define INT_SIMD          __m512i
#       define RI_SIMD           INT_SIMD

#   endif


    /* The "SIMD_BLOCK" value can be played with.  It has no effect on the
     * accuracy, but affects the performance.  Too low or too high values can
     * decrease the computation speed. */
#   define SIMD_BLOCK 4


#   define MUL_R(x, y)         PF(mul)((x), (y))
#   define DIV_R(x, y)         PF(div)((x), (y))
#   define ADD_R(x, y)         PF(add)((x), (y))
#   define SUB_R(x, y)         PF(sub)((x), (y))


#   define SET1_R(x)           PF(set1)((x))
#   define LOAD_R(x)           PF(load)((x))
#   define LOADU_R(x)          PF(loadu)((x))
#   define STORE_R(ptr, x)     PF(store)((ptr), (x))
#   define STOREU_R(ptr, x)    PF(storeu)((ptr), (x))


#   if HAVE_AVX || HAVE_AVX2
#       define BLEND_R(x, y, mask)   PF(blendv)((x), (y), (mask))
#       define MOVEMASK(x)           PF(movemask)((x))
#   elif HAVE_AVX512F
#       define BLEND_R(x, y, mask)   PF(mask_blend)((mask), (x), (y))
#       define MOVEMASK(x)           (x)
#   endif


#   define AND_R(x, y)          PF(and)((x), (y))
#   define OR_R(x, y)           PF(or)((x), (y))
#   if HAVE_AVX || HAVE_AVX2
#       define EQ_R(x, y)       PF(cmp)((x), (y), _CMP_EQ_OQ)
#       define GE_R(x, y)       PF(cmp)((x), (y), _CMP_GE_OQ)
#       define LT_R(x, y)       PF(cmp)((x), (y), _CMP_LT_OQ)
#   elif HAVE_AVX512F
#       define EQ_R(x, y)       PF_MASK(cmp)((x), (y), _CMP_EQ_OQ)
#       define GE_R(x, y)       PF_MASK(cmp)((x), (y), _CMP_GE_OQ)
#       define LT_R(x, y)       PF_MASK(cmp)((x), (y), _CMP_LT_OQ)
#   endif


#   if HAVE_AVX

#       define SET1_RI                SET1_R
#       if CHARM_FLOAT
#           define LOAD_RI(x)         PF(set)((x)[7], (x)[6], (x)[5], \
                                              (x)[4], (x)[3], (x)[2], \
                                              (x)[1], (x)[0])
#       else
#           define LOAD_RI(x)         PF(set)((x)[3], (x)[2], \
                                              (x)[1], (x)[0])
#       endif
#       define EQ_RI                  EQ_R
#       define GT_RI(x, y)            PF(cmp)((x), (y), _CMP_GT_OQ)
#       define ADD_RI                 ADD_R
#       define SUB_RI                 SUB_R
#       define CAST_RI2R(x)           (x)
#       define CAST_R2RI(x)           (x)
#       define OR_MASK                OR_R
#       define SET_ZERO_R             PF(setzero)()
#       define SET_ZERO_I             P(setzero_si256)()
#       define SET_ZERO_RI            SET_ZERO_R
#       define ANDNOT_MASK(x)         EQ_R((x), SET_ZERO_R)
#       define BLEND_RI               BLEND_R

#   elif HAVE_AVX2

#       define SET1_RI(x)             PINT2(set1)((x))
#       if CHARM_FLOAT
#           define LOAD_RI(x)         PINT2(set)((x)[7], (x)[6], \
                                                 (x)[5], (x)[4], \
                                                 (x)[3], (x)[2], \
                                                 (x)[1], (x)[0])
#           define CAST_R2RI(x)       P(castps_si256)((x))
#       else
#           define LOAD_RI(x)         PINT2(set)((x)[3], (x)[2], \
                                                 (x)[1], (x)[0])
#           define CAST_R2RI(x)       P(castpd_si256)((x))
#       endif
#       define EQ_RI(x, y)            PINT(cmpeq)((x), (y))
#       define GT_RI(x, y)            PINT(cmpgt)((x), (y))
#       define ADD_RI(x, y)           PINT(add)((x), (y))
#       define SUB_RI(x, y)           PINT(sub)((x), (y))
#       define CAST_RI2R(x)           PF(castsi256)((x))
#       define AND_RI(x, y)           P(and_si256)((x), (y))
#       define OR_MASK(x, y)          P(or_si256)((x), (y))
#       define SET_ZERO_R             PF(setzero)()
#       define SET_ZERO_I             P(setzero_si256)()
#       define SET_ZERO_RI            SET_ZERO_I
#       define ANDNOT_MASK(x)         EQ_RI((x), SET_ZERO_RI)
#       define BLEND_RI(x, y, mask)   CAST_R2RI(BLEND_R(CAST_RI2R(x), \
                                                        CAST_RI2R(y), \
                                                        CAST_RI2R(mask)))
#   elif HAVE_AVX512F

#       define SET1_RI(x)             PINT(set1)((x))
#       if CHARM_FLOAT
#           define LOAD_RI(x)         P(set_epi32)((x)[15], (x)[14], \
                                                   (x)[13], (x)[12], \
                                                   (x)[11], (x)[10], \
                                                   (x)[9],  (x)[8],  \
                                                   (x)[7],  (x)[6],  \
                                                   (x)[5],  (x)[4],  \
                                                   (x)[3],  (x)[2],  \
                                                   (x)[1],  (x)[0])
#           define OR_MASK(x, y)      _kor_mask16((x), (y))
#       else
#           define LOAD_RI(x)         PINT(set)((x)[7], (x)[6], \
                                                (x)[5], (x)[4], \
                                                (x)[3], (x)[2], \
                                                (x)[1], (x)[0])
#           define OR_MASK(x, y)      _kor_mask8((x), (y))
#       endif
#       define EQ_RI(x, y)            PI_MASK(cmpeq)((x), (y))
#       define GT_RI(x, y)            PI_MASK(cmpgt)((x), (y))
#       define ADD_RI(x, y)           PINT(add)((x), (y))
#       define SUB_RI(x, y)           PINT(sub)((x), (y))
#       define CAST_RI2R(x)           (x)
#       define CAST_R2RI(x)           (x)
#       define SET_ZERO_R             PF(setzero)()
#       define SET_ZERO_I             P(setzero_si512)()
#       define SET_ZERO_RI            SET_ZERO_I
#       define ANDNOT_MASK(x)         _knot_mask8((x))
#       define BLEND_RI(x, y, mask)   PINT(mask_blend)((mask), (x), (y))

#   endif


    /* Absolute value of a SIMD vector */
    /* ..................................................................... */
    /* Compute the absolute value of all elements of a SIMD vector using the
     * "ABS_R" macro.  Before using "ABS_R" in a code, the "ABS_R_INIT" macro
     * must be called, ideally only once outside any loop. */
#   ifdef CHARM_FLOAT
#       define NONSIGNBITS 0x7FFFFFFF
#   else
#       define NONSIGNBITS 0x7FFFFFFFFFFFFFFFLL
#   endif
#   if HAVE_AVX || HAVE_AVX2
#       define ABS_R_INIT   REAL_SIMD NONSIGNBITS_R = \
                                        PF(castsi256)(PINT2(set1)(NONSIGNBITS))
#   elif HAVE_AVX512F
#       define ABS_R_INIT   REAL_SIMD NONSIGNBITS_R = \
                                        PF(castsi512)(PINT2(set1)(NONSIGNBITS))
#   endif
#   define ABS_R(x)     PF(and)(NONSIGNBITS_R, (x))
    /* ..................................................................... */



    /* Change the sign of a SIMD vector */
    /* ..................................................................... */
    /* Similarly as with the "ABS_R" macro, also "NEG_R_INIT" must be called
     * before calling "NEG_R", ideally only once outside any loop. */
#   ifdef CHARM_FLOAT
#       define SIGNBIT 0x80000000
#   else
#       define SIGNBIT 0x8000000000000000LL
#   endif
#   if HAVE_AVX || HAVE_AVX2
#       define NEG_R_INIT   REAL_SIMD SIGNBIT_R = \
                                        PF(castsi256)(PINT2(set1)(SIGNBIT))
#   elif HAVE_AVX512F
#       define NEG_R_INIT   REAL_SIMD SIGNBIT_R = \
                                        PF(castsi512)(PINT2(set1)(SIGNBIT))
#   endif
#   define NEG_R(x)     PF(xor)((x), SIGNBIT_R)
    /* ..................................................................... */


    /* Get the smallest multiple of "x" that is equal to or larger
     * than "multiple". */
#   define SIMD_MULTIPLE(x, multiple) ((((x) + \
                                         (multiple) - 1) / (multiple)) * \
                                       (multiple))


    /* Sum all elements of a SIMD vector. */
#   if HAVE_AVX || HAVE_AVX2
        inline static REAL SUM_R(REAL_SIMD x)
        {
            x = ADD_R(PF(permute2f128)(x, x, 0x01), x);
#       if CHARM_FLOAT
            x = PF(hadd)(x, x);
            x = PF(hadd)(x, x);
            return P(cvtss_f32)(x);
#       else
            x = ADD_R(PF(permute)(x, 1), x);
            return P(cvtsd_f64)(x);
#       endif
        }
#   elif HAVE_AVX512F
#       define SUM_R(x)    PF(reduce_add)((x))
#   endif


#else


#   define SIMD_SIZE     1
    /* The "SIMD_BLOCK" value can be played with.  It has no effect on the
     * accuracy, but affects the performance.  Too low or too high values can
     * decrease the computation speed. */
#   define SIMD_BLOCK    8
#   define SIMD_TRUE     1
#   define SIMD_MEMALIGN 0
#   define REAL_SIMD     REAL
#   define INT_SIMD      int
#   define RI_SIMD       int


#   define MUL_R(x, y)    ((x) * (y))
#   define DIV_R(x, y)    ((x) / (y))
#   define ADD_R(x, y)    ((x) + (y))
#   define SUB_R(x, y)    ((x) - (y))


#   define SET1_R(x)            (x)
#   define SET1_RI(x)           (x)
#   define SET_ZERO_R           (PREC(0.0))
#   define SET_ZERO_I           (0)
#   define SET_ZERO_RI          (0)
#   define STORE_R(ptr, x)      (*(ptr) = (x))
#   define STOREU_R             STORE_R
#   define LOAD_R(x)            (*(x))
#   define LOADU_R(x)           (*(x))
#   define SUM_R(x)             (x)


#   define EQ_R(x, y)           ((x) == (y))
#   define LT_R(x, y)           ((x) < (y))
#   define MOVEMASK(x)          (x)


#   define MASK_SIMD            int
#   define MASK2_SIMD           MASK_SIMD


    /* Absolute value.  The "FABS" macro is defined in "../prec.h". */
#   define ABS_R_INIT
#   define NONSIGNBITS
#   define ABS_R                FABS


    /* Change the sign. */
#   define NEG_R_INIT
#   define SIGNBIT
#   define NEG_R(x)             (-(x))


#   define SIMD_MULTIPLE(x, multiple) (x)


#endif


#endif
