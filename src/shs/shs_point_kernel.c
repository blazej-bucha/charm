/* ------------------------------------------------------------------------- */
/* This file implements the sums over harmonic degrees that appear in spherical
 * harmonic synthesis.  Its main job is to implement routines to compute the
 * lumped coefficients, depending on the order of the radial and latitudinal
 * derivatives.  In general, computed are the following quantities:
 *
 *      a_m = \sum_{n = m}^{nmax} ampl(n) * (R / r)^{n + 1 + i}
 *            * \bar{P}_{nm}(\sin\varphi) * \bar{C}_{nm},
 *
 *      b_m = \sum_{n = m}^{nmax} ampl(n) * (R / r)^{n + 1 + i}
 *            * \bar{P}_{nm}(\sin\varphi) * \bar{S}_{nm},
 *
 * where
 *
 *      "i" is the order of the potential derivative ("0" for the gravitational
 *      potential, "1" for the gravitational vector, "2" for the gravitational
 *      tensor),
 *
 *      "ampl(n)" is amplitude factor, which depends on the order of the
 *      potential derivative ("1" if "i == 0", "n + 1" if "i == 1", "(n + 1)
 *      * (n + 2)" if "i == 2").
 *
 * These equations are computed if "KERNEL_GRAD == 0".  If "KERNEL_GRAD == 1",
 * computed are on-the-fly the three groups of lumped coefficients for all
 * three derivatives with respect to "r", "phi" and "lambda".  If "KERNEL_GRAD
 * == 2", six groups of lumped coefficients are computed accordingly.
 *
 * The outputs are stored in the "lc" variable.  If the synthesis is being
 * computed for a symmetric grid, the "a_m" and "b_m" quantities for the
 * latitude "-\varphi" will be stored in "lc->a2" and "lc->b2", respectively.
 * The same rule holds for other members of "lc".
 *
 * This source file heavily relies on MACROS, so it is rather difficult to read
 * and debug.  However, MACROs make it possible to avoid (1) unnecessary "if"
 * statements in performance-critical for loops and (2) overheads from calling
 * external routines.  These factors could otherwise seriously degrade
 * performance.
 *
 * One should not try to compile directly this code.  Instead, this source file
 * is included several times in "shs_point_kernels.c", each time with
 * a different value of the symbolic constants "DR" (order of the derivative
 * with respect to the radius), "DLAT" (order of the derivative with respect to
 * the latitude) and "DLON" (order of the longitudinal derivative).  As
 * a result, we get kernels for each combination of "DR", "DLAT", "DLON", etc.
 *
 * WARNING: Be extremely careful when doing any modifications of this file. */
/* ------------------------------------------------------------------------- */






/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "../shc/shc_block_struct.h"
#include "../shc/shc_block_get_idx.h"
#include "../simd/simd.h"
#include "../leg/leg_func_xnum.h"
#include "../leg/leg_func_use_xnum.h"
#include "../misc/misc_sd_calloc.h"
#include "../glob/glob_get_shs_block_lat_multiplier.h"
#include "shs_check_symm_simd.h"
#include "shs_point_gradn.h"
#include "shs_lc_struct.h"
#include "shs_point_kernels.h"
/* ------------------------------------------------------------------------- */






/* This source file is included in "shs_point_kernels.c" and should not be
 * compiled directly.  To avoid this, we compile it only "#if COMPILE_KERNELS"
 * is true. */
#if COMPILE_KERNELS


/* ------------------------------------------------------------------------- */
/* "DR" represents the order of the potential derivative in the radial
 * direction. */
#ifndef DR
#   error "DR undefined."
#endif
#if (DR != 0) && (DR != 1) && (DR != 2) && (DR != GRAD_1) && (DR != GRAD_2)
#   error "DR must be an integer -2, -1, ..., 2."
#endif


/* "DLAT" represents the order of the potential derivative in the latitudinal
 * direction. */
#ifndef DLAT
#   error "DLAT undefined."
#endif
#if (DLAT != 0) && (DLAT != 1) && (DLAT != 2) && (DLAT != GRAD_1) && \
    (DLAT != GRAD_2)
#   error "DLAT must be an integer -2, -1, ..., 2."
#endif


/* "DLON" represents the order of the potential derivative in the longitudinal
 * direction. */
#ifndef DLON
#   error "DLON undefined."
#endif
#if (DLON != 0) && (DLON != 1) && (DLON != 2) && (DLON != GRAD_1) && \
    (DLON != GRAD_2)
#   error "DLON must be an integer -2, -1, ..., 2."
#endif


/* If "KERNEL_GRAD" is "1", compiled will be the kernel for the full
 * first-order gradient.  If "2", compiled will be the kernel for the full
 * second-order gradient.  If "0", compiled will be a kernel for a single
 * zero-, first- or second-order derivative based on the values of "DR", "DLAT"
 * and "DLON". */
#undef KERNEL_GRAD


#if (DR == GRAD_1) && (DLAT == GRAD_1) && (DLON == GRAD_1)
    /* We are asked to compile the kernel for the full first-order gradient, so
     * set "KERNEL_GRAD" to "1" to indicate this choice later. */
#   undef KERNEL_GRAD
#   define KERNEL_GRAD 1
#   undef DR  /* "DR" is now "-1", so make it "1" as needed below. */
#   define DR 1
#   undef DLAT  /* "DLAT" is now "-1", so make it "1" as needed below. */
#   define DLAT 1
#   undef DLON  /* "DLON" is now "-1", so make it "1" as needed below. */
#   define DLON 1
#elif (DR == GRAD_2) && (DLAT == GRAD_2) && (DLON == GRAD_2)
    /* The same but for the second-order gradient. */
#   undef KERNEL_GRAD
#   define KERNEL_GRAD 2
#   undef DR  /* "DR" is now "-2", so make it "2" as needed below. */
#   define DR 2
#   undef DLAT  /* "DLAT" is now "-2", so make it "2" as needed below. */
#   define DLAT 2
#   undef DLON  /* "DLON" is now "-2", so make it "2" as needed below. */
#   define DLON 2
#else
#   define KERNEL_GRAD 0
#endif


/* gcc 15.2.1 miscompiles kernels from this source file if SIMD is disabled and
 * "-O3 -ffast-math" (or "-march=native", etc.) are enabled.  It seems that gcc
 * incorrectly reorders some statements that must not be reordered, thereby
 * producing grossly invalid results.  "MEMORY_BARRIER" ensures that gcc will
 * never reorder the statements found before and after the macro.
 * "MEMORY_BARRIER" is non-empty only if compiling with gcc, SIMD is disabled
 * and gcc version 4.4 or newer is used (required by the assembly call). */
#undef MEMORY_BARRIER
#define MEMORY_BARRIER
#if defined(__GNUC__) && !defined(SIMD)
#   if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 4)
#       undef MEMORY_BARRIER
#       define MEMORY_BARRIER asm volatile("" : : : "memory");
#   endif
#endif
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
/* The "AMPL" macro defines the "1", "(n + 1)" or "(n + 1) * (n + 2)" amplitude
 * factors that depend on the order of the radial derivative of the
 * gravitational potential. */
#undef AMPL
#if KERNEL_GRAD == 0
#   if DR == 0
        /* "1" */
#       define AMPL(n)                                                        \
            ;
#   elif DR == 1
        /* (n + 1) */
#       define AMPL(n)                                                        \
               ampl = SET1_R(((n) + 1));
#   elif (DR == 2)
        /* (n + 1) * (n + 2) */
#       define AMPL(n)                                                        \
               ampl = SET1_R(((n) + 1) * ((n) + 2));
#   endif
#elif KERNEL_GRAD == 1
        /* (n + 1) */
#   define AMPL(n)                                                            \
           ampl1 = SET1_R(((n) + 1));
#elif KERNEL_GRAD == 2
    /* ampl1 = (n + 1), ampl2 = (n + 2)*/
#   define AMPL(n)                                                            \
           ampl1 = SET1_R(((n) + 1));                                         \
           ampl2 = SET1_R(((n) + 2));
#endif


/* Get the offset for "rpows" and "rpows2" to get the suitable power "(R
 * / r)^(n + 1 + p)" for a given degree "n" and the "p"th derivative of the
 * "potential" */
#undef GET_RPOWS_IDX
#define GET_RPOWS_IDX(n) ((n + 1 + dorder) * BLOCK_S)


/* When implementing the latitudinal symmetry property of Legendre functions,
 * we end up with equations where the sums for "a2" and "b2" have alternating
 * signs with respect to degree "n".  To allow an easy implementation, that is,
 * without writing similar code multiple times, we define the "SIGN1" and
 * "SIGN2" macros.  Based on "SIGN1" or "SIGN2", we either add or subtract the
 * "n"th-degree contribution to "a2" and "b2".  This approach is especially
 * useful, given that the signs for the first-order derivatives of Legendre
 * functions are switched with respect to the Legendre functions themselves.
 * So here, we specify whether the "n"th-degree contribution is added or
 * subracted when applying the symmetry of Legendre functions and of their
 * derivatives. */
#undef SIGN1
#undef SIGN2
#if KERNEL_GRAD == 0
#   if DLAT == 0
#      define SIGN1 SUB_R
#      define SIGN2 ADD_R
#   elif DLAT == 1
#      define SIGN1 ADD_R
#      define SIGN2 SUB_R
#   elif DLAT == 2
#      define SIGN1 SUB_R
#      define SIGN2 ADD_R
#   endif
#else
#   define SIGN1 SUB_R
#   define SIGN2 ADD_R
#endif


/* We stored the Legendre functions in variables "pnm0", "pnm1" and "pnm2"
 * ("P_{n - 2, m}", "P_{n - 1, m}", "P_{n, m}", respectively).  The first-order
 * derivatives of "p_{n, m}" are stored in "dpnm0", "dpnm1", "dpnm2" and the
 * second-order derivatives in "ddpnm0", "ddpnm1", "ddpnm2".  It is handy to
 * introduce the "DIFF" macro that adds the "d" or "dd" prefixes (or nothing)
 * before "pnm0", depending on the order of the derivative for which we compile
 * the kernel.  This significantly simplifies the code. */
#undef DIFF
#if DLAT == 0
#   define DIFF(x)   x
#elif DLAT == 1
#   define DIFF(x)   CAT(d, x)
#elif DLAT == 2
#   define DIFF(x)   CAT(dd, x)
#endif


/* This macro computes the first-order derivatives of Legendre functions for
 * all latitudes within a "BLOCK_S" if "DLAT > 0".  Otherwise, it does
 * nothing. */
#undef DPNM_RECURRENCE_BLOCK
#if DLAT > 0
#   define DPNM_RECURRENCE_BLOCK                                              \
        for (l = 0; l < BLOCK_S; l++)                                         \
        {                                                                     \
            DPNM_RECURRENCE((dpnm2[l]), (x[l]), (pnm2[l]), (tu[l]),           \
                            (u_rec[l]), (ns), (enms));                        \
        }
#else
#   define DPNM_RECURRENCE_BLOCK
#endif


/* This macro computes the second-order derivatives of Legendre functions for
 * all latitudes within a "BLOCK_S" if "DLAT > 1".  Otherwise, it does
 * nothing. */
#undef DDPNM_RECURRENCE_BLOCK
#if DLAT > 1
#   define DDPNM_RECURRENCE_BLOCK                                             \
        for (l = 0; l < BLOCK_S; l++)                                         \
        {                                                                     \
            DDPNM_RECURRENCE((ddpnm2[l]), (dpnm2[l]), (pnm2[l]),              \
                             (tu[l]), (u2_rec[l]), (m2s), (nn1s));            \
        }
#else
#   define DDPNM_RECURRENCE_BLOCK
#endif


/* Lumped coefficients for "grad0" */
#undef GRAD0_LC
#undef GRAD0_LC_R1
#if KERNEL_GRAD == 0

#   define GRAD0_LC(PM_R, ab, i, cs)                                          \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT(ab, i)[l] = PM_R(lc->CAT(ab, i)[l],                    \
                                        MUL_R(CAT2(rpows, i, _m)[l],          \
                                              CAT2(leg_, cs, nm)[l]));        \
           }

#   define GRAD0_LC_R1(PM_R, ab, i, cs)                                       \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT(ab, i)[l] = PM_R(lc->CAT(ab, i)[l],                    \
                                        CAT2(leg_, cs, nm)[l]);               \
           }
#else
#   define GRAD0_LC(PM_R, ab, i, cs) ;
#   define GRAD0_LC_R1(PM_R, ab, i, cs) ;
#endif


/* Lumped coefficients for "grad1" */
#undef GRAD1_LC
#undef GRAD1_LC_R1
#if KERNEL_GRAD > 0

#   define GRAD1_LC(PM1_R, PM2_R, ab, i, cs)                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT(ab, i)[l] = PM1_R(lc->CAT(ab, i)[l],                   \
                                         MUL_R(CAT2(rpows, i, _m)[l],         \
                                               CAT2(leg_, cs, nm)[l]));       \
           }                                                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, r, i)[l] = PM1_R(lc->CAT2(ab, r, i)[l],           \
                                             MUL_R(CAT2(rpows, i, _m)[l],     \
                                                   CAT2(leg_, cs, nm_r)[l])); \
           }                                                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, p, i)[l] = PM2_R(lc->CAT2(ab, p, i)[l],           \
                                             MUL_R(CAT2(rpows, i, _m)[l],     \
                                                   CAT2(leg_, cs, nm_p)[l])); \
           }

#   define GRAD1_LC_R1(PM1_R, PM2_R, ab, i, cs)                               \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT(ab, i)[l] = PM1_R(lc->CAT(ab, i)[l],                   \
                                         CAT2(leg_, cs, nm)[l]);              \
           }                                                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, r, i)[l] = PM1_R(lc->CAT2(ab, r, i)[l],           \
                                             CAT2(leg_, cs, nm_r)[l]);        \
           }                                                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, p, i)[l] = PM2_R(lc->CAT2(ab, p, i)[l],           \
                                             CAT2(leg_, cs, nm_p)[l]);        \
           }
#else
#   define GRAD1_LC(PM1_R, PM2_R, ab, i, cs) ;
#   define GRAD1_LC_R1(PM1_R, PM2_R, ab, i, cs) ;
#endif


/* Lumped coefficients for "grad2" */
#undef GRAD2_LC
#undef GRAD2_LC_R1
#if KERNEL_GRAD > 1

#   define GRAD2_LC(PM1_R, PM2_R, ab, i, cs)                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, rr, i)[l] = PM1_R(lc->CAT2(ab, rr, i)[l],         \
                                              MUL_R(CAT2(rpows, i, _m)[l],    \
                                                  CAT2(leg_, cs, nm_rr)[l])); \
           }                                                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, rp, i)[l] = PM2_R(lc->CAT2(ab, rp, i)[l],         \
                                              MUL_R(CAT2(rpows, i, _m)[l],    \
                                                  CAT2(leg_, cs, nm_rp)[l])); \
           }                                                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, pp, i)[l] = PM1_R(lc->CAT2(ab, pp, i)[l],         \
                                              MUL_R(CAT2(rpows, i, _m)[l],    \
                                                  CAT2(leg_, cs, nm_pp)[l])); \
           }

#   define GRAD2_LC_R1(PM1_R, PM2_R, ab, i, cs)                               \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, rr, i)[l] = PM1_R(lc->CAT2(ab, rr, i)[l],         \
                                              CAT2(leg_, cs, nm_rr)[l]);      \
           }                                                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, rp, i)[l] = PM2_R(lc->CAT2(ab, rp, i)[l],         \
                                              CAT2(leg_, cs, nm_rp)[l]);      \
           }                                                                  \
           for (l = 0; l < BLOCK_S; l++)                                      \
           {                                                                  \
               lc->CAT2(ab, pp, i)[l] = PM1_R(lc->CAT2(ab, pp, i)[l],         \
                                              CAT2(leg_, cs, nm_pp)[l]);      \
           }
#else
#   define GRAD2_LC(PM1_R, PM2_R, ab, i, cs) ;
#   define GRAD2_LC_R1(PM1_R, PM2_R, ab, i, cs) ;
#endif


/* Products of Legendre functions (or their derivatives) and the spherical
 * harmonic coefficients for "grad0" */
#undef GRAD0_LEGCS
#if KERNEL_GRAD == 0
#   if DR == 0
#       define GRAD0_LEGCS(pnm, cs)                                           \
              for (l = 0; l < BLOCK_S; l++)                                   \
              {                                                               \
                  CAT2(leg_, cs, nm)[l] = MUL_R(DIFF(pnm[l]), CAT(cs, nm));   \
              }
#   else
#       define GRAD0_LEGCS(pnm, cs)                                           \
              for (l = 0; l < BLOCK_S; l++)                                   \
              {                                                               \
                  CAT2(leg_, cs, nm)[l] = MUL_R(DIFF(pnm[l]), CAT(cs, nm));   \
                  CAT2(leg_, cs, nm)[l] = MUL_R(ampl, CAT2(leg_, cs, nm)[l]); \
              }
#   endif
#else
#   define GRAD0_LEGCS(pnm, cs) ;
#endif


/* Products of Legendre functions (or their derivatives) and the spherical
 * harmonic coefficients for "grad1" */
#undef GRAD1_LEGCS
#if KERNEL_GRAD > 0
#   define GRAD1_LEGCS(pnm, dpnm, cs)                                         \
        for (l = 0; l < BLOCK_S; l++)                                         \
        {                                                                     \
            CAT2(leg_, cs, nm)[l] = MUL_R(pnm[l], CAT(cs, nm));               \
        }                                                                     \
        for (l = 0; l < BLOCK_S; l++)                                         \
        {                                                                     \
            CAT2(leg_, cs, nm_r)[l] = MUL_R(ampl1, CAT2(leg_, cs, nm)[l]);    \
        }                                                                     \
        for (l = 0; l < BLOCK_S; l++)                                         \
        {                                                                     \
            CAT2(leg_, cs, nm_p)[l] = MUL_R(dpnm[l], CAT(cs, nm));            \
        }
#else
#   define GRAD1_LEGCS(pnm, dpnm, cs) ;
#endif


/* Products of Legendre functions (or their derivatives) and the spherical
 * harmonic coefficients for "grad2" */
#undef GRAD2_LEGCS
#if KERNEL_GRAD > 1
#   define GRAD2_LEGCS(pnm, cs)                                               \
        for (l = 0; l < BLOCK_S; l++)                                         \
        {                                                                     \
            CAT2(leg_, cs, nm_rr)[l] = MUL_R(ampl2, CAT2(leg_, cs, nm_r)[l]); \
        }                                                                     \
        for (l = 0; l < BLOCK_S; l++)                                         \
        {                                                                     \
            CAT2(leg_, cs, nm_rp)[l] = MUL_R(ampl1, CAT2(leg_, cs, nm_p)[l]); \
        }                                                                     \
        for (l = 0; l < BLOCK_S; l++)                                         \
        {                                                                     \
            CAT2(leg_, cs, nm_pp)[l] = MUL_R(pnm[l], CAT(cs, nm));            \
        }
#else
#   define GRAD2_LEGCS(pnm, cs) ;
#endif


/* A single macro for all products of Legendre functions (or their derivatives)
 * and the spherical harmonic coefficients */
#undef LFCS
#define LFCS                                                                  \
    GRAD0_LEGCS(pnm2, c);                                                     \
    GRAD0_LEGCS(pnm2, s);                                                     \
    GRAD1_LEGCS(pnm2, dpnm2, c);                                              \
    GRAD1_LEGCS(pnm2, dpnm2, s);                                              \
    GRAD2_LEGCS(ddpnm2, c);                                                   \
    GRAD2_LEGCS(ddpnm2, s);


/* Similarly, a single macro for lumped coefficients, all grads and "(R / r)
 * > 1" */
#undef LCAB
#define LCAB(n, PM1_R, PM2_R)                                                 \
    rpows_m = rpows + GET_RPOWS_IDX(n);                                       \
    GRAD0_LC(ADD_R, a, , c);                                                  \
    GRAD0_LC(ADD_R, b, , s);                                                  \
    GRAD1_LC(ADD_R, ADD_R, a, , c);                                           \
    GRAD1_LC(ADD_R, ADD_R, b, , s);                                           \
    GRAD2_LC(ADD_R, ADD_R, a, , c);                                           \
    GRAD2_LC(ADD_R, ADD_R, b, , s);                                           \
                                                                              \
                                                                              \
    if (symm)                                                                 \
    {                                                                         \
        rpows2_m = rpows2 + GET_RPOWS_IDX(n);                                 \
        GRAD0_LC(PM1_R, a, 2, c);                                             \
        GRAD0_LC(PM1_R, b, 2, s);                                             \
        GRAD1_LC(PM1_R, PM2_R, a, 2, c);                                      \
        GRAD1_LC(PM1_R, PM2_R, b, 2, s);                                      \
        GRAD2_LC(PM1_R, PM2_R, a, 2, c);                                      \
        GRAD2_LC(PM1_R, PM2_R, b, 2, s);                                      \
    }




/* Similarly, a single macro for lumped coefficients, all grads and "(R / r) ==
 * 1" */
#undef LCAB_R1
#define LCAB_R1(PM1_R, PM2_R)                                                 \
    GRAD0_LC_R1(ADD_R, a, , c);                                               \
    GRAD0_LC_R1(ADD_R, b, , s);                                               \
    GRAD1_LC_R1(ADD_R, ADD_R, a, , c);                                        \
    GRAD1_LC_R1(ADD_R, ADD_R, b, , s);                                        \
    GRAD2_LC_R1(ADD_R, ADD_R, a, , c);                                        \
    GRAD2_LC_R1(ADD_R, ADD_R, b, , s);                                        \
                                                                              \
                                                                              \
    if (symm)                                                                 \
    {                                                                         \
        GRAD0_LC_R1(PM1_R, a, 2, c);                                          \
        GRAD0_LC_R1(PM1_R, b, 2, s);                                          \
        GRAD1_LC_R1(PM1_R, PM2_R, a, 2, c);                                   \
        GRAD1_LC_R1(PM1_R, PM2_R, b, 2, s);                                   \
        GRAD2_LC_R1(PM1_R, PM2_R, a, 2, c);                                   \
        GRAD2_LC_R1(PM1_R, PM2_R, b, 2, s);                                   \
    }


#undef ENMS
#if DLAT > 0
#   define ENMS(n)                                                            \
        enms = SET1_R(enm[(n)]);                                              \
        ns   = SET1_R((REAL)(n));
#else
#   define ENMS(n)
#endif


#undef NN1S
#if DLAT > 1
#   define NN1S(n)                                                            \
        nn1s = SET1_R((REAL)((n) * ((n) + 1)));
#else
#   define NN1S(n)
#endif


/* A single macro for Legendre functions and their derivatives */
#undef LEG_CS
#define LEG_CS(n)                                                             \
    anms = SET1_R(anm[(n)]);                                                  \
    bnms = SET1_R(bnm[(n)]);                                                  \
    ENMS(n);                                                                  \
    NN1S(n);                                                                  \
    cnm  = SET1_R(shcs_block->c[idx]);                                        \
    snm  = SET1_R(shcs_block->s[idx]);                                        \
    AMPL((n));                                                                \
                                                                              \
                                                                              \
    for (l = 0; l < BLOCK_S; l++)                                             \
    {                                                                         \
        PNM_RECURRENCE(x[l], y[l], pnm2[l], t[l], anms, bnms);                \
    }                                                                         \
    DPNM_RECURRENCE_BLOCK;                                                    \
    DDPNM_RECURRENCE_BLOCK;                                                   \
                                                                              \
                                                                              \
    for (l = 0; l < BLOCK_S; l++)                                             \
    {                                                                         \
        RECURRENCE_NEXT_ITER(y[l], x[l], pnm2[l]);                            \
    }                                                                         \
                                                                              \
                                                                              \
    LFCS


/* The "LOOP_ITER" macro can be used for any values of "R / r".
 *
 * The "LOOP_ITER_R1" macro can *only* be used if all "(R / r) == 1.0".
 *
 * "LOOP_ITER_R1" is faster than "LOOP_ITER", so the "is_ratio_one" value is
 * used in the main code to check whether "LOOP_ITER_R1" can safely be used or
 * not.
 *
 * Except for "LOOP_ITER_R1", the entire kernel assumes that "(R / r) != 1", so
 * that the slower code variant is used elsewhere.  This include code to
 * compute Legendre polynomials, sectorial Legendre functions, and tesseral
 * Legendre functions (the last one with X-numbers).  This is because the
 * computing time of these parts of the code are rather negligible and using
 * a separate code for "(R / r) == 1.0" would make the code uselessly
 * lengthy. */
#define LOOP_ITER(n, PM1_R, PM2_R)                                            \
    LEG_CS((n))                                                               \
    LCAB((n), PM1_R, PM2_R)                                                   \
    idx++;


#define LOOP_ITER_R1(n, PM1_R, PM2_R)                                         \
    LEG_CS((n))                                                               \
    LCAB_R1(PM1_R, PM2_R)                                                     \
    idx++;
/* ------------------------------------------------------------------------- */






/* ------------------------------------------------------------------------- */
#ifndef COMPILE_STATIC_FUNCTIONS
#define COMPILE_STATIC_FUNCTIONS


/* Swaps pointers "x" and "y". */
static inline void swap(REAL_SIMD **x,
                        REAL_SIMD **y,
                        REAL_SIMD *tmp)
{
    tmp = *x;
    *x  = *y;
    *y  = tmp;
    return;
}






/* Change sign of all elements of "x". */
static inline void change_sign(REAL_SIMD *x,
                               size_t n)
{
    for (size_t i = 0; i < n; i++)
        x[i] = NEG_R(x[i]);


    return;
}






/* Multiply all elements of "x" and "y" by "m_factor" (used to compute the
 * longitudinal derivatives). */
static inline void apply_m_factor(REAL_SIMD m_factor,
                                  REAL_SIMD *x,
                                  REAL_SIMD *y,
                                  size_t n)
{
    for (size_t i = 0; i < n; i++)
        x[i] = MUL_R(m_factor, x[i]);


    for (size_t i = 0; i < n; i++)
        y[i] = MUL_R(m_factor, y[i]);


    return;
}


/* First-order longitudinal derivative. */
static inline void dlon1(REAL_SIMD **a,
                         REAL_SIMD **b,
                         REAL_SIMD **a2,
                         REAL_SIMD **b2,
                         unsigned long m,
                         _Bool symm,
                         size_t n)
{
    REAL_SIMD *tmp = NULL;
    swap(a, b, tmp);
    change_sign(*b, n);
    REAL_SIMD m_factor = SET1_R((REAL)m);
    apply_m_factor(m_factor, *a, *b, n);


    if (symm)
    {
        swap(a2, b2, tmp);
        change_sign(*b2, n);
        apply_m_factor(m_factor, *a2, *b2, n);
    }


    return;
}






/* Second-order longitudinal derivative. */
static inline void dlon2(REAL_SIMD *a,
                         REAL_SIMD *b,
                         REAL_SIMD *a2,
                         REAL_SIMD *b2,
                         unsigned long m,
                         _Bool symm,
                         size_t n)
{
    change_sign(a, n);
    change_sign(b, n);
    REAL_SIMD m_factor = SET1_R((REAL)(m * m));
    apply_m_factor(m_factor, a, b, n);


    if (symm)
    {
        change_sign(a2, n);
        change_sign(b2, n);
        apply_m_factor(m_factor, a2, b2, n);
    }


    return;
}


#endif
/* ------------------------------------------------------------------------- */






/* Returns summations over harmonic degrees, "lc", for synthesis with point
 * values */
#if KERNEL_GRAD == 0
#   if   (DR == 0) && (DLAT == 0) && (DLON == 0)
void CHARM(shs_point_kernel_dr0_dlat0_dlon0)
#   elif (DR == 1) && (DLAT == 0) && (DLON == 0)
void CHARM(shs_point_kernel_dr1_dlat0_dlon0)
#   elif (DR == 2) && (DLAT == 0) && (DLON == 0)
void CHARM(shs_point_kernel_dr2_dlat0_dlon0)
#   elif (DR == 0) && (DLAT == 1) && (DLON == 0)
void CHARM(shs_point_kernel_dr0_dlat1_dlon0)
#   elif (DR == 0) && (DLAT == 2) && (DLON == 0)
void CHARM(shs_point_kernel_dr0_dlat2_dlon0)
#   elif (DR == 0) && (DLAT == 0) && (DLON == 1)
void CHARM(shs_point_kernel_dr0_dlat0_dlon1)
#   elif (DR == 0) && (DLAT == 0) && (DLON == 2)
void CHARM(shs_point_kernel_dr0_dlat0_dlon2)
#   elif (DR == 1) && (DLAT == 1) && (DLON == 0)
void CHARM(shs_point_kernel_dr1_dlat1_dlon0)
#   elif (DR == 1) && (DLAT == 0) && (DLON == 1)
void CHARM(shs_point_kernel_dr1_dlat0_dlon1)
#   elif (DR == 0) && (DLAT == 1) && (DLON == 1)
void CHARM(shs_point_kernel_dr0_dlat1_dlon1)
#   else
#       error "Wrong combination of DR, DLAT and DLON."
#   endif
#elif KERNEL_GRAD == 1
void CHARM(shs_point_kernel_grad1)
#elif KERNEL_GRAD == 2
void CHARM(shs_point_kernel_grad2)
#else
#   error "Wrong value of KERNEL_GRAD."
#endif
                            (unsigned long nmax,
                             unsigned long m,
                             const CHARM(shc_block) *shcs_block,
                             _Bool is_ratio_one,
                             const REAL *anm,
                             const REAL *bnm,
                             const REAL *enm,
                             const REAL_SIMD *t,
                             const REAL_SIMD *u,
                             const REAL *ps,
                             const INT *ips,
                             REAL_SIMD *rpows,
                             REAL_SIMD *rpows2,
                             const REAL_SIMD *symm_simd,
                             unsigned dorder,
                             CHARM(lc) *lc)
{
    /* --------------------------------------------------------------------- */
#if HAVE_MPI
    const size_t BLOCK_S = CHARM(glob_get_shs_block_lat_multiplier)();
#else
#   define BLOCK_S SIMD_BLOCK_S
#endif


    MISC_SD_CALLOC_REAL_SIMD_INIT(x);
    MISC_SD_CALLOC_REAL_SIMD_INIT(y);
    MISC_SD_CALLOC_REAL_SIMD_INIT(z);
    MISC_SD_CALLOC_RI_SIMD_INIT(ix);
    MISC_SD_CALLOC_RI_SIMD_INIT(iy);
    MISC_SD_CALLOC_RI_SIMD_INIT(iz);
    MISC_SD_CALLOC_RI_SIMD_INIT(ixy);
    MISC_SD_CALLOC_REAL_SIMD_INIT(pnm0);
    MISC_SD_CALLOC_REAL_SIMD_INIT(pnm1);
    MISC_SD_CALLOC_REAL_SIMD_INIT(pnm2);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_cnm);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_snm);
    MISC_SD_CALLOC__BOOL_INIT(ds);
#if DLAT > 0
    MISC_SD_CALLOC_REAL_SIMD_INIT(tu);
    MISC_SD_CALLOC_REAL_SIMD_INIT(dpnm0);
    MISC_SD_CALLOC_REAL_SIMD_INIT(dpnm1);
    MISC_SD_CALLOC_REAL_SIMD_INIT(dpnm2);
#endif
#if DLAT > 1
    MISC_SD_CALLOC_REAL_SIMD_INIT(ddpnm0);
    MISC_SD_CALLOC_REAL_SIMD_INIT(ddpnm1);
    MISC_SD_CALLOC_REAL_SIMD_INIT(ddpnm2);
#endif
#if (DLAT > 0) || (DLON > 0)
    MISC_SD_CALLOC_REAL_SIMD_INIT(u_rec);
#endif
#if (DLAT > 1) || (DLON > 1)
    MISC_SD_CALLOC_REAL_SIMD_INIT(u2_rec);
#endif
#if KERNEL_GRAD > 0
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_cnm_r);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_snm_r);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_cnm_p);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_snm_p);
#endif
#if KERNEL_GRAD > 1
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_cnm_rr);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_snm_rr);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_cnm_rp);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_snm_rp);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_cnm_pp);
    MISC_SD_CALLOC_REAL_SIMD_INIT(leg_snm_pp);
#endif


    MISC_SD_CALLOC_REAL_SIMD_E(x, BLOCK_S, SIMD_BLOCK_S, lc->error, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(y, BLOCK_S, SIMD_BLOCK_S, lc->error, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(z, BLOCK_S, SIMD_BLOCK_S, lc->error, BARRIER_1);
    MISC_SD_CALLOC_RI_SIMD_E(ix, BLOCK_S, SIMD_BLOCK_S, lc->error, BARRIER_1);
    MISC_SD_CALLOC_RI_SIMD_E(iy, BLOCK_S, SIMD_BLOCK_S, lc->error, BARRIER_1);
    MISC_SD_CALLOC_RI_SIMD_E(iz, BLOCK_S, SIMD_BLOCK_S, lc->error, BARRIER_1);
    MISC_SD_CALLOC_RI_SIMD_E(ixy, BLOCK_S, SIMD_BLOCK_S, lc->error, BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(pnm0, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(pnm1, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(pnm2, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_cnm, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_snm, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC__BOOL_E(ds, BLOCK_S, SIMD_BLOCK_S, lc->error, BARRIER_1);
#if DLAT > 0
    MISC_SD_CALLOC_REAL_SIMD_E(tu, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(dpnm0, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(dpnm1, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(dpnm2, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
#endif
#if DLAT > 1
    MISC_SD_CALLOC_REAL_SIMD_E(ddpnm0, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(ddpnm1, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(ddpnm2, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
#endif
#if (DLAT > 0) || (DLON > 0)
    MISC_SD_CALLOC_REAL_SIMD_E(u_rec, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
#endif
#if (DLAT > 1) || (DLON > 1)
    MISC_SD_CALLOC_REAL_SIMD_E(u2_rec, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
#endif
#if KERNEL_GRAD > 0
    MISC_SD_CALLOC_REAL_SIMD_E(leg_cnm_r, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_snm_r, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_cnm_p, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_snm_p, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
#endif
#if KERNEL_GRAD > 1
    MISC_SD_CALLOC_REAL_SIMD_E(leg_cnm_rr, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_snm_rr, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_cnm_rp, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_snm_rp, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_cnm_pp, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
    MISC_SD_CALLOC_REAL_SIMD_E(leg_snm_pp, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_1);
#endif


#if HAVE_MPI
BARRIER_1:
    if (lc->error)
        goto FAILURE_1;
#endif
    /* --------------------------------------------------------------------- */


    REAL_SIMD w;
#ifdef SIMD
    const RI_SIMD    zero_ri = SET_ZERO_RI;
    const RI_SIMD    one_ri  = SET1_RI(1);
    const RI_SIMD    mone_ri = SET1_RI(-1);
    const REAL_SIMD  zero_r  = SET_ZERO_R;
    const REAL_SIMD  BIG_r   = SET1_R(BIG);
    const REAL_SIMD  BIGI_r  = SET1_R(BIGI);
    const REAL_SIMD  BIGS_r  = SET1_R(BIGS);
    const REAL_SIMD  BIGSI_r = SET1_R(BIGSI);
    REAL_SIMD  tmp1_r,  tmp2_r;
    MASK_SIMD  mask1, mask2;
    MASK2_SIMD mask3;
#endif
    const REAL_SIMD one = SET1_R(PREC(1.0));


    size_t l;
#if (DLAT > 0) || (DLON > 0)
    for (l = 0; l < BLOCK_S; l++)
        u_rec[l] = DIV_R(one, u[l]);  /* 1 / cos(lat) */
#endif


#if DLAT > 0
    for (l = 0; l < BLOCK_S; l++)
        tu[l] = MUL_R(t[l], u_rec[l]);  /* sin(lat) / cos(lat) */
#endif
#if (DLAT > 1) || (DLON > 1)
    for (l = 0; l < BLOCK_S; l++)
        u2_rec[l] = MUL_R(u_rec[l], u_rec[l]);  /* 1 / cos(lat)^2 */
#endif
    const REAL_SIMD ROOT3_r = SET1_R(ROOT3);


    REAL_SIMD cnm, snm;


    _Bool npm_even; /* True if "n + m" is even */
    _Bool symm = 0;
    for (l = 0; l < BLOCK_S; l++)
        symm = symm || CHARM(shs_check_symm_simd)(symm_simd[l]);


    /* Get the index of "Cmm" and "Smm" in "shcs_block->c" and "shcs_block->s"
     * */
    unsigned long idx = CHARM(shc_block_get_idx)(shcs_block, m);


    REAL_SIMD *rpows_m;
    REAL_SIMD *rpows2_m;


    REAL_SIMD anms, bnms;
#if DLAT > 0
    REAL_SIMD enms;
    REAL_SIMD ns;  /* "n" */
#endif
#if DLAT > 1
    REAL_SIMD nn1s;  /* "n * (n + 1)" */
    REAL_SIMD m2s = SET1_R((REAL)(m * m));  /* "m * m" */
#endif
#if KERNEL_GRAD == 0
#   if DR > 0
    REAL_SIMD ampl;
#   endif
#else
    REAL_SIMD ampl1;
#   if KERNEL_GRAD > 1
    REAL_SIMD ampl2;
#   endif
#endif


    /* Reset "a2" and "b2".  Required. */
    for (l = 0; l < BLOCK_S; l++)
        lc->a2[l] = lc->b2[l] = SET_ZERO_R;
#if KERNEL_GRAD > 0
    for (l = 0; l < BLOCK_S; l++)
        lc->ar2[l] = lc->br2[l] = lc->ap2[l] = lc->bp2[l] = SET_ZERO_R;
#endif
#if KERNEL_GRAD > 1
    for (l = 0; l < BLOCK_S; l++)
        lc->arr2[l] = lc->brr2[l] = lc->arp2[l] = lc->brp2[l] =
            lc->app2[l] = lc->bpp2[l] = SET_ZERO_R;
#endif


    /* Summations over harmonic degree "n" */
    if (m == 0)
    {
        /* Zonal harmonics */
        /* ----------------------------------------------------- */
        /* P00 */
        cnm = SET1_R(shcs_block->c[idx]);
        idx++;
        AMPL(0);


        for (l = 0; l < BLOCK_S; l++)
            pnm0[l] = one;
#if DLAT > 0
        for (l = 0; l < BLOCK_S; l++)
            dpnm0[l] = SET_ZERO_R;
#endif
#if DLAT > 1
        for (l = 0; l < BLOCK_S; l++)
            ddpnm0[l] = SET_ZERO_R;
#endif


#if KERNEL_GRAD == 0
        GRAD0_LEGCS(pnm0, c);
#endif
#if KERNEL_GRAD > 0
        GRAD1_LEGCS(pnm0, dpnm0, c);
#endif
#if KERNEL_GRAD > 1
        GRAD2_LEGCS(ddpnm0, c);
#endif
        rpows_m = rpows + GET_RPOWS_IDX(0);
        for (l = 0; l < BLOCK_S; l++)
            lc->a[l] = MUL_R(rpows_m[l], leg_cnm[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->b[l] = SET_ZERO_R;
#if KERNEL_GRAD > 0
        for (l = 0; l < BLOCK_S; l++)
            lc->ar[l] = MUL_R(rpows_m[l], leg_cnm_r[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->br[l] = SET_ZERO_R;
        for (l = 0; l < BLOCK_S; l++)
            lc->ap[l] = MUL_R(rpows_m[l], leg_cnm_p[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->bp[l] = SET_ZERO_R;
#endif
#if KERNEL_GRAD > 1
        for (l = 0; l < BLOCK_S; l++)
            lc->arr[l] = MUL_R(rpows_m[l], leg_cnm_rr[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->arp[l] = MUL_R(rpows_m[l], leg_cnm_rp[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->app[l] = MUL_R(rpows_m[l], leg_cnm_pp[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->brr[l] = lc->brp[l] = lc->bpp[l] = SET_ZERO_R;
#endif


        if (symm)
        {
            rpows2_m = rpows2 + GET_RPOWS_IDX(0);
#if KERNEL_GRAD == 0
            GRAD0_LC(SIGN2, a, 2, c);
            for (l = 0; l < BLOCK_S; l++)
                lc->b2[l] = SET_ZERO_R;
#endif
#if KERNEL_GRAD > 0
            GRAD1_LC(ADD_R, ADD_R, a, 2, c);
            for (l = 0; l < BLOCK_S; l++)
                lc->b2[l] = lc->br2[l] = lc->bp2[l] = SET_ZERO_R;
#endif
#if KERNEL_GRAD > 1
            GRAD2_LC(ADD_R, SUB_R, a, 2, c);
            for (l = 0; l < BLOCK_S; l++)
                lc->brr2[l] = lc->brp2[l] = lc->bpp2[l] = SET_ZERO_R;
#endif
        }


        /* P10 */
        if (nmax >= 1)
        {
            cnm = SET1_R(shcs_block->c[idx]);
            idx++;
            AMPL(1);


            for (l = 0; l < BLOCK_S; l++)
                pnm1[l] = MUL_R(ROOT3_r, t[l]);
#if DLAT > 0
            for (l = 0; l < BLOCK_S; l++)
                dpnm1[l] = MUL_R(ROOT3_r, u[l]);
#endif
#if DLAT > 1
            for (l = 0; l < BLOCK_S; l++)
                ddpnm1[l] = NEG_R(pnm1[l]);
#endif


#if KERNEL_GRAD == 0
            GRAD0_LEGCS(pnm1, c);
#endif
#if KERNEL_GRAD > 0
            GRAD1_LEGCS(pnm1, dpnm1, c);
#endif
#if KERNEL_GRAD > 1
            GRAD2_LEGCS(ddpnm1, c);
#endif
            rpows_m = rpows + GET_RPOWS_IDX(1);
#if KERNEL_GRAD == 0
            GRAD0_LC(ADD_R, a, , c);
#endif
#if KERNEL_GRAD > 0
            GRAD1_LC(ADD_R, ADD_R, a, , c);
#endif
#if KERNEL_GRAD > 1
            GRAD2_LC(ADD_R, ADD_R, a, , c);
#endif


            if (symm)
            {
                rpows2_m = rpows2 + GET_RPOWS_IDX(1);
#if KERNEL_GRAD == 0
                GRAD0_LC(SIGN1, a, 2, c);
#endif
#if KERNEL_GRAD > 0
                GRAD1_LC(SUB_R, ADD_R, a, 2, c);
#endif
#if KERNEL_GRAD > 1
                GRAD2_LC(SUB_R, ADD_R, a, 2, c);
#endif
            }
        }


        /* P20, P30, ..., Pnmax,0 */
        if (nmax >= 2)
        {
            /* Is "n + m" even?  Since we start the loop with "n = 2" and "m
             * = 0", then the parity of the first "n + m" is always even.
             * Then, it changes with every loop iteration. */
            npm_even = 1;


            for (unsigned long n = 2; n <= nmax; n++, npm_even = !npm_even)
            {
                anms = SET1_R(anm[n]);
                bnms = SET1_R(bnm[n]);
#if DLAT > 0
                enms = SET1_R(enm[n]);
                ns   = SET1_R((REAL)n);
#endif
#if DLAT > 1
                nn1s = SET1_R((REAL)(n * (n + 1)));
#endif
                cnm  = SET1_R(shcs_block->c[idx]);
                idx++;
                AMPL(n);


                for (l = 0; l < BLOCK_S; l++)
                    PNM_RECURRENCE(pnm1[l], pnm0[l], pnm2[l], t[l], anms,
                                   bnms);
#if DLAT > 0
                for (l = 0; l < BLOCK_S; l++)
                    DPNM_RECURRENCE(dpnm2[l], pnm1[l], pnm2[l], tu[l],
                                    u_rec[l], ns, enms);
#endif
#if DLAT > 1
                for (l = 0; l < BLOCK_S; l++)
                    DDPNM_RECURRENCE(ddpnm2[l], dpnm2[l], pnm2[l], tu[l],
                                     u2_rec[l], m2s, nn1s);
#endif


#if KERNEL_GRAD == 0
                GRAD0_LEGCS(pnm2, c);
#endif
#if KERNEL_GRAD > 0
                GRAD1_LEGCS(pnm2, dpnm2, c);
#endif
#if KERNEL_GRAD > 1
                GRAD2_LEGCS(ddpnm2, c);
#endif


                rpows_m = rpows + GET_RPOWS_IDX(n);
#if KERNEL_GRAD == 0
                GRAD0_LC(ADD_R, a, , c);
#endif
#if KERNEL_GRAD > 0
                GRAD1_LC(ADD_R, ADD_R, a, , c);
#endif
#if KERNEL_GRAD > 1
                GRAD2_LC(ADD_R, ADD_R, a, , c);
#endif


                for (l = 0; l < BLOCK_S; l++)
                    RECURRENCE_NEXT_ITER(pnm0[l], pnm1[l], pnm2[l]);


                if (symm)
                {
                    rpows2_m = rpows2 + GET_RPOWS_IDX(n);
                    if (npm_even)
                    {
#if KERNEL_GRAD == 0
                        GRAD0_LC(SIGN2, a, 2, c);
#endif
#if KERNEL_GRAD > 0
                        GRAD1_LC(ADD_R, SUB_R, a, 2, c);
#endif
#if KERNEL_GRAD > 1
                        GRAD2_LC(ADD_R, SUB_R, a, 2, c);
#endif
                    }
                    else
                    {
#if KERNEL_GRAD == 0
                        GRAD0_LC(SIGN1, a, 2, c);
#endif
#if KERNEL_GRAD > 0
                        GRAD1_LC(SUB_R, ADD_R, a, 2, c);
#endif
#if KERNEL_GRAD > 1
                        GRAD2_LC(SUB_R, ADD_R, a, 2, c);
#endif
                    }
                }
            }
        }
        /* ----------------------------------------------------- */
    }
    else /* Non-zonal harmonics */
    {

        /* Sectorial harmonics */
        /* ----------------------------------------------------- */
        cnm = SET1_R(shcs_block->c[idx]);
        snm = SET1_R(shcs_block->s[idx]);
        idx++;
#if DLAT > 0
        enms = SET1_R(enm[m]);
        ns   = SET1_R((REAL)m);
#endif
#if DLAT > 1
        nn1s = SET1_R(m * (m + 1));
#endif
        AMPL(m);


        for (l = 0; l < BLOCK_S; l++)
        {
#ifdef SIMD
            PNM_SECTORIAL_XNUM_SIMD(x[l], ix[l],
                                    ps[(SIMD_SIZE * nmax) * l + (m - 1) *
                                       SIMD_SIZE],
                                    ips[(SIMD_SIZE * nmax) * l + (m - 1) *
                                        SIMD_SIZE],
                                    pnm0[l], BIG_r, zero_r, zero_ri, mone_ri,
                                    mask1, mask2, SECTORIALS);
#else
            PNM_SECTORIAL_XNUM(x[l], ix[l],
                               ps[(SIMD_SIZE * nmax) * l + (m - 1) *
                                  SIMD_SIZE],
                               ips[(SIMD_SIZE * nmax) * l + (m - 1) *
                                   SIMD_SIZE], pnm0[l]);
#endif
        }
#if DLAT > 0
        for (l = 0; l < BLOCK_S; l++)
            dpnm0[l] = NEG_R(MUL_R(MUL_R(MUL_R(ns, t[l]), u_rec[l]), pnm0[l]));
#endif
#if DLAT > 1
        for (l = 0; l < BLOCK_S; l++)
            DDPNM_RECURRENCE(ddpnm0[l], dpnm0[l], pnm0[l], tu[l], u2_rec[l],
                             m2s, nn1s);
#endif


#if KERNEL_GRAD == 0
        GRAD0_LEGCS(pnm0, c);
        GRAD0_LEGCS(pnm0, s);
#endif
#if KERNEL_GRAD > 0
        GRAD1_LEGCS(pnm0, dpnm0, c);
        GRAD1_LEGCS(pnm0, dpnm0, s);
#endif
#if KERNEL_GRAD > 1
        GRAD2_LEGCS(ddpnm0, c);
        GRAD2_LEGCS(ddpnm0, s);
#endif


        rpows_m = rpows + GET_RPOWS_IDX(m);
        for (l = 0; l < BLOCK_S; l++)
            lc->a[l] = MUL_R(rpows_m[l], leg_cnm[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->b[l] = MUL_R(rpows_m[l], leg_snm[l]);
#if KERNEL_GRAD > 0
        for (l = 0; l < BLOCK_S; l++)
            lc->ar[l] = MUL_R(rpows_m[l], leg_cnm_r[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->ap[l] = MUL_R(rpows_m[l], leg_cnm_p[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->br[l] = MUL_R(rpows_m[l], leg_snm_r[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->bp[l] = MUL_R(rpows_m[l], leg_snm_p[l]);
#endif
#if KERNEL_GRAD > 1
        for (l = 0; l < BLOCK_S; l++)
            lc->arr[l] = MUL_R(rpows_m[l], leg_cnm_rr[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->arp[l] = MUL_R(rpows_m[l], leg_cnm_rp[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->app[l] = MUL_R(rpows_m[l], leg_cnm_pp[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->brr[l] = MUL_R(rpows_m[l], leg_snm_rr[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->brp[l] = MUL_R(rpows_m[l], leg_snm_rp[l]);
        for (l = 0; l < BLOCK_S; l++)
            lc->bpp[l] = MUL_R(rpows_m[l], leg_snm_pp[l]);
#endif


        if (symm)
        {
            rpows2_m = rpows2 + GET_RPOWS_IDX(m);
#if KERNEL_GRAD == 0
            GRAD0_LC(SIGN2, a, 2, c);
            GRAD0_LC(SIGN2, b, 2, s);
#endif
#if KERNEL_GRAD > 0
            GRAD1_LC(ADD_R, SUB_R, a, 2, c);
            GRAD1_LC(ADD_R, SUB_R, b, 2, s);
#endif
#if KERNEL_GRAD > 1
            GRAD2_LC(ADD_R, SUB_R, a, 2, c);
            GRAD2_LC(ADD_R, SUB_R, b, 2, s);
#endif
        }
        /* ----------------------------------------------------- */


        /* Tesseral harmonics */
        /* ----------------------------------------------------- */
        if (m < nmax)
        {
            anms = SET1_R(anm[m + 1]);
            bnms = SET1_R(bnm[m + 1]);
#if DLAT > 0
            enms = SET1_R(enm[m + 1]);
            ns   = SET1_R((REAL)(m + 1));
#endif
#if DLAT > 1
            nn1s = SET1_R((m + 1) * (m + 2));
#endif
            cnm  = SET1_R(shcs_block->c[idx]);
            snm  = SET1_R(shcs_block->s[idx]);
            idx++;
            AMPL(m + 1);


            for (l = 0; l < BLOCK_S; l++)
            {
#ifdef SIMD
                PNM_SEMISECTORIAL_XNUM_SIMD(x[l], y[l], ix[l], iy[l], w, t[l],
                                            anms,
                                            pnm1[l], mask1, mask2, mask3,
                                            zero_r, zero_ri, mone_ri, BIG_r,
                                            BIGS_r,  BIGI_r, SEMISECTORIALS);
#else
                PNM_SEMISECTORIAL_XNUM(x[l], y[l], ix[l], iy[l], w, t[l],
                                       anms, pnm1[l]);
#endif
            }
#if DLAT > 0
            for (l = 0; l < BLOCK_S; l++)
                DPNM_RECURRENCE(dpnm1[l], pnm0[l], pnm1[l], tu[l], u_rec[l],
                                ns, enms);
#endif
#if DLAT > 1
            for (l = 0; l < BLOCK_S; l++)
                DDPNM_RECURRENCE(ddpnm1[l], dpnm1[l], pnm1[l], tu[l],
                                 u2_rec[l], m2s, nn1s);
#endif


#if KERNEL_GRAD == 0
            GRAD0_LEGCS(pnm1, c);
            GRAD0_LEGCS(pnm1, s);
#endif
#if KERNEL_GRAD > 0
            GRAD1_LEGCS(pnm1, dpnm1, c);
            GRAD1_LEGCS(pnm1, dpnm1, s);
#endif
#if KERNEL_GRAD > 1
            GRAD2_LEGCS(ddpnm1, c);
            GRAD2_LEGCS(ddpnm1, s);
#endif


            rpows_m = rpows + GET_RPOWS_IDX(m + 1);
#if KERNEL_GRAD == 0
            GRAD0_LC(ADD_R, a, , c);
            GRAD0_LC(ADD_R, b, , s);
#endif
#if KERNEL_GRAD > 0
            GRAD1_LC(ADD_R, ADD_R, a, , c);
            GRAD1_LC(ADD_R, ADD_R, b, , s);
#endif
#if KERNEL_GRAD > 1
            GRAD2_LC(ADD_R, ADD_R, a, , c);
            GRAD2_LC(ADD_R, ADD_R, b, , s);
#endif


            if (symm)
            {
                rpows2_m = rpows2 + GET_RPOWS_IDX(m + 1);
#if KERNEL_GRAD == 0
                GRAD0_LC(SIGN1, a, 2, c);
                GRAD0_LC(SIGN1, b, 2, s);
#endif
#if KERNEL_GRAD > 0
                GRAD1_LC(SUB_R, ADD_R, a, 2, c);
                GRAD1_LC(SUB_R, ADD_R, b, 2, s);
#endif
#if KERNEL_GRAD > 1
                GRAD2_LC(SUB_R, ADD_R, a, 2, c);
                GRAD2_LC(SUB_R, ADD_R, b, 2, s);
#endif
            }


            /* Loop over degrees */
            /* ------------------------------------------------- */
            for (l = 0; l < BLOCK_S; l++)
                ds[l] = 0;


            /* Is "n + m" even?  Since we start the loop with "n = m + 2", then
             * the parity of the first "m + 2 + m" is always even.  Then, it
             * changes with every loop iteration. */
            npm_even = 1;


            unsigned long n;
            for (n = (m + 2);
                 CHARM(leg_func_use_xnum(ds, BLOCK_S)) && n <= nmax;
                 n++, npm_even = !npm_even)
            {
                anms = SET1_R(anm[n]);
                bnms = SET1_R(bnm[n]);
#if DLAT > 0
                enms = SET1_R(enm[n]);
                ns   = SET1_R((REAL)n);
#endif
#if DLAT > 1
                nn1s = SET1_R((REAL)(n * (n + 1)));
#endif
                cnm  = SET1_R(shcs_block->c[idx]);
                snm  = SET1_R(shcs_block->s[idx]);
                idx++;
                AMPL(n);


                /* Compute tesseral Legendre function */
                for (l = 0; l < BLOCK_S; l++)
                {
#ifdef SIMD
                    PNM_TESSERAL_XNUM_SIMD(x[l], y[l], z[l],
                                           ix[l], iy[l], iz[l], ixy[l],
                                           w, t[l], anms, bnms, pnm2[l],
                                           tmp1_r, tmp2_r, mask1, mask2,
                                           mask3, zero_r, zero_ri, one_ri,
                                           BIG_r, BIGI_r, BIGS_r, BIGSI_r,
                                           TESSERALS1, TESSERALS2, ds[l]);
#else
                    PNM_TESSERAL_XNUM(x[l], y[l], z[l], ix[l], iy[l], iz[l],
                                      ixy[l], w, t[l], anms, bnms, pnm2[l],
                                      ds[l]);
#endif
                }
#if DLAT > 0
                for (l = 0; l < BLOCK_S; l++)
                    DPNM_RECURRENCE(dpnm2[l], pnm1[l], pnm2[l], tu[l],
                                    u_rec[l], ns, enms);
                for (l = 0; l < BLOCK_S; l++)
                    pnm1[l]  = pnm2[l];
#endif
#if DLAT > 1
                for (l = 0; l < BLOCK_S; l++)
                    DDPNM_RECURRENCE(ddpnm2[l], dpnm2[l], pnm2[l], tu[l],
                                     u2_rec[l], m2s, nn1s);
#endif


#if KERNEL_GRAD == 0
                GRAD0_LEGCS(pnm2, c);
                GRAD0_LEGCS(pnm2, s);
#endif
#if KERNEL_GRAD > 0
                GRAD1_LEGCS(pnm2, dpnm2, c);
                GRAD1_LEGCS(pnm2, dpnm2, s);
#endif
#if KERNEL_GRAD > 1
                GRAD2_LEGCS(ddpnm2, c);
                GRAD2_LEGCS(ddpnm2, s);
#endif


                rpows_m = rpows + GET_RPOWS_IDX(n);
#if KERNEL_GRAD == 0
                GRAD0_LC(ADD_R, a, , c);
                GRAD0_LC(ADD_R, b, , s);
#endif
#if KERNEL_GRAD > 0
                GRAD1_LC(ADD_R, ADD_R, a, , c);
                GRAD1_LC(ADD_R, ADD_R, b, , s);
#endif
#if KERNEL_GRAD > 1
                GRAD2_LC(ADD_R, ADD_R, a, , c);
                GRAD2_LC(ADD_R, ADD_R, b, , s);
#endif


                if (symm)
                {
                    rpows2_m = rpows2 + GET_RPOWS_IDX(n);
                    if (npm_even)
                    {
#if KERNEL_GRAD == 0
                        GRAD0_LC(SIGN2, a, 2, c);
                        GRAD0_LC(SIGN2, b, 2, s);
#endif
#if KERNEL_GRAD > 0
                        GRAD1_LC(ADD_R, SUB_R, a, 2, c);
                        GRAD1_LC(ADD_R, SUB_R, b, 2, s);
#endif
#if KERNEL_GRAD > 1
                        GRAD2_LC(ADD_R, SUB_R, a, 2, c);
                        GRAD2_LC(ADD_R, SUB_R, b, 2, s);
#endif
                    }
                    else
                    {
#if KERNEL_GRAD == 0
                        GRAD0_LC(SIGN1, a, 2, c);
                        GRAD0_LC(SIGN1, b, 2, s);
#endif
#if KERNEL_GRAD > 0
                        GRAD1_LC(SUB_R, ADD_R, a, 2, c);
                        GRAD1_LC(SUB_R, ADD_R, b, 2, s);
#endif
#if KERNEL_GRAD > 1
                        GRAD2_LC(SUB_R, ADD_R, a, 2, c);
                        GRAD2_LC(SUB_R, ADD_R, b, 2, s);
#endif
                    }
                }
            }


            if (n > nmax)
                goto DR_DERIVATIVE;


            /* From now on, "F"-numbers can be used instead of the "X"-numbers
             * to gain some speed */


            /* We want to unroll the loop that follows.  To do that, we need to
             * make sure that the loop starts with an even value of "n + m".
             * So if "n + m" is odd, we need to hard code one iteration. */
            if (!npm_even)
            {
                LOOP_ITER(n, SIGN1, SIGN2);
                MEMORY_BARRIER;
                n++;
            }


            if (is_ratio_one)
            {
                for (; (n + 1) <= nmax; n += 2)
                {
                    LOOP_ITER_R1(n,     SIGN2, SIGN1);
                    MEMORY_BARRIER;
                    LOOP_ITER_R1(n + 1, SIGN1, SIGN2);
                }
            }
            else
            {
                for (; (n + 1) <= nmax; n += 2)
                {
                    LOOP_ITER(n,     SIGN2, SIGN1);
                    MEMORY_BARRIER;
                    LOOP_ITER(n + 1, SIGN1, SIGN2);
                }
            }


            if (n > nmax)
                goto DR_DERIVATIVE;

            MEMORY_BARRIER;
            LOOP_ITER(n, SIGN2, SIGN1);
            /* ------------------------------------------------- */


        } /* End of computation of tesseral harmonics */
        /* ----------------------------------------------------- */


    } /* End of the summation over harmonic degree "n" */
    /* --------------------------------------------------------- */






DR_DERIVATIVE:
    ;  /* Required to avoid a gcc warning: "a label can only be part of
        * a statement and a declaration is not a statement". */


    /* For the lumped coefficients of the first-order derivative, we need to
     * change their sign right. */
    /* --------------------------------------------------------- */
#if (KERNEL_GRAD == 0) && (DR == 1)
    change_sign(lc->a, BLOCK_S);
    change_sign(lc->b, BLOCK_S);


    if (symm)
    {
        change_sign(lc->a2, BLOCK_S);
        change_sign(lc->b2, BLOCK_S);
    }
#endif
#if KERNEL_GRAD > 0
    change_sign(lc->ar, BLOCK_S);
    change_sign(lc->br, BLOCK_S);


    if (symm)
    {
        change_sign(lc->ar2, BLOCK_S);
        change_sign(lc->br2, BLOCK_S);
    }
#endif
#if KERNEL_GRAD > 1
    change_sign(lc->arp, BLOCK_S);
    change_sign(lc->brp, BLOCK_S);


    if (symm)
    {
        change_sign(lc->arp2, BLOCK_S);
        change_sign(lc->brp2, BLOCK_S);
    }
#endif
    /* --------------------------------------------------------- */






    /* --------------------------------------------------------- */
#if KERNEL_GRAD == 0  /* Longitudinal derivatives if computing a single
                       * parameter only */

    /* Apply the "1 / cos(lat)" or "1 / cos^2(lat)" terms if needed */
    /* ..................................................................... */
#   undef CLAT_TERM
#   if DLON == 1

#       define CLAT_TERM(x)                                                   \
            for (l = 0; l < BLOCK_S; l++)                                     \
                lc->x[l] = MUL_R(u_rec[l], lc->x[l]);

#   elif DLON == 2

#       define CLAT_TERM(x)                                                   \
            for (l = 0; l < BLOCK_S; l++)                                     \
                lc->x[l] = MUL_R(u2_rec[l], lc->x[l]);

#   else

#       define CLAT_TERM(x)

#   endif


    CLAT_TERM(a);
    CLAT_TERM(b);
    if (symm)
    {
        CLAT_TERM(a2);
        CLAT_TERM(b2);
    }
    /* ..................................................................... */


    /* Longitudinal derivatives */
    /* ..................................................................... */
#   if DLON == 1

    dlon1(&lc->a, &lc->b, &lc->a2, &lc->b2, m, symm, BLOCK_S);

#   elif DLON == 2

    dlon2(lc->a, lc->b, lc->a2, lc->b2, m, symm, BLOCK_S);

#   endif
    /* ..................................................................... */


#elif KERNEL_GRAD == 1  /* The full first-order gradient in LNOF */

    /* Longitudinal derivative */
    dlon1(&lc->a, &lc->b, &lc->a2, &lc->b2, m, symm, BLOCK_S);

    /* In LNOF, we need negative first-order longitudinal derivative, so change
     * the signs in "lc->a" and "lc->b". */
    change_sign(lc->a, BLOCK_S);
    change_sign(lc->b, BLOCK_S);
    if (symm)
    {
        change_sign(lc->a2, BLOCK_S);
        change_sign(lc->b2, BLOCK_S);
    }


    /* Now divide the lumped coefficients by the "cos(lat)" term. */
#undef Y
#define Y(x)                                                                  \
    for (l = 0; l < BLOCK_S; l++)                                             \
        lc->x[l] = MUL_R(u_rec[l], lc->x[l]);


    Y(a);
    Y(b);
    if (symm)
    {
        Y(a2);
        Y(b2);
    }

#elif KERNEL_GRAD == 2  /* The full second-order gradient in LNOF */

    /* xx */
    /* ..................................................................... */
#undef XX
#define XX(ab, i)                                                             \
    for (l = 0; l < BLOCK_S; l++)                                             \
        lc->CAT2(ab, pp, i)[l] = ADD_R(lc->CAT2(ab, r, i)[l],                 \
                                       lc->CAT2(ab, pp, i)[l]);


    XX(a, );
    XX(b, );


    if (symm)
    {
        XX(a, 2);
        XX(b, 2);
    }
    /* ..................................................................... */


    /* yy */
    /* ..................................................................... */
    /* We need a hard copy of "lc->a", "lc->b", "lc->a2" and "lc->b2", because
     * "dlon2" overwrites these with the second order longitudinal derivatives
     * and later we need also the first order derivatives, which need to be
     * derived from the original variables. */
    MISC_SD_CALLOC_REAL_SIMD_INIT(a_);
    MISC_SD_CALLOC_REAL_SIMD_INIT(b_);
    MISC_SD_CALLOC_REAL_SIMD_INIT(a_2);
    MISC_SD_CALLOC_REAL_SIMD_INIT(b_2);


    MISC_SD_CALLOC_REAL_SIMD_E(a_, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_2);
    MISC_SD_CALLOC_REAL_SIMD_E(b_, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_2);
    MISC_SD_CALLOC_REAL_SIMD_E(a_2, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_2);
    MISC_SD_CALLOC_REAL_SIMD_E(b_2, BLOCK_S, SIMD_BLOCK_S, lc->error,
                               BARRIER_2);


#if HAVE_MPI
BARRIER_2:
    if (lc->error)
        goto FAILURE_2;
#endif


    for (l = 0; l < BLOCK_S; l++)
        a_[l] = lc->a[l];
    for (l = 0; l < BLOCK_S; l++)
        b_[l] = lc->b[l];
    if (symm)
    {
        for (l = 0; l < BLOCK_S; l++)
            a_2[l] = lc->a2[l];
        for (l = 0; l < BLOCK_S; l++)
            b_2[l] = lc->b2[l];
    }
    REAL_SIMD *aptr_  = &a_[0];
    REAL_SIMD *bptr_  = &b_[0];
    REAL_SIMD *aptr_2 = &a_2[0];
    REAL_SIMD *bptr_2 = &b_2[0];


    dlon2(lc->a, lc->b, lc->a2, lc->b2, m, symm, BLOCK_S);


    /* If "symm", we need to change the sign before the "Vp" term, hence the
     * "PM" input parameter. */
#undef YY
#define YY(ab, PM, i)                                                         \
    for (l = 0; l < BLOCK_S; l++)                                             \
        lc->CAT(ab, i)[l] = ADD_R(PM(MUL_R(u2_rec[l],                         \
                                           lc->CAT(ab, i)[l]),                \
                                     MUL_R(tu[l], lc->CAT2(ab, p, i)[l])),    \
                                  lc->CAT2(ab, r, i)[l]);
    YY(a, SUB_R, );
    YY(b, SUB_R, );


    if (symm)
    {
        YY(a, ADD_R, 2);
        YY(b, ADD_R, 2);
    }
    /* ..................................................................... */


    /* xz */
    /* ..................................................................... */
#undef XZ
#define XZ(ab, i)                                                             \
    for (l = 0; l < BLOCK_S; l++)                                             \
        lc->CAT2(ab, rp, i)[l] = SUB_R(lc->CAT2(ab, rp, i)[l],                \
                                       lc->CAT2(ab, p, i)[l]);


    XZ(a, );
    XZ(b, );


    if (symm)
    {
        XZ(a, 2);
        XZ(b, 2);
    }
    /* ..................................................................... */


    /* xy */
    /* ..................................................................... */
    /* Vlp */
    dlon1(&lc->ap, &lc->bp, &lc->ap2, &lc->bp2, m, symm, BLOCK_S);
    /* Vl */
    dlon1(&aptr_, &bptr_, &aptr_2, &bptr_2, m, symm, BLOCK_S);


#undef XY
#define XY(ab, PM, i)                                                         \
    for (l = 0; l < BLOCK_S; l++)                                             \
        lc->CAT2(ab, p, i)[l] = PM(NEG_R(MUL_R(u_rec[l],                      \
                                               lc->CAT2(ab, p, i)[l])),       \
                                   MUL_R(MUL_R(tu[l], u_rec[l]),              \
                                         CAT2(ab, ptr_, i)[l]));


    XY(a, SUB_R, );
    XY(b, SUB_R, );


    if (symm)
    {
        XY(a, ADD_R, 2);
        XY(b, ADD_R, 2);
    }
    /* ..................................................................... */


    /* yz */
    /* ..................................................................... */
    /* Vlp */
    dlon1(&lc->ar, &lc->br, &lc->ar2, &lc->br2, m, symm, BLOCK_S);


#undef YZ
#define YZ(ab, i)                                                             \
    for (l = 0; l < BLOCK_S; l++)                                             \
        lc->CAT2(ab, r, i)[l] = MUL_R(u_rec[l],                               \
                                      SUB_R(CAT2(ab, ptr_, i)[l],             \
                                            lc->CAT2(ab, r, i)[l]));


    YZ(a, );
    YZ(b, );


    if (symm)
    {
        YZ(a, 2);
        YZ(b, 2);
    }
    /* ..................................................................... */


#if HAVE_MPI
FAILURE_2:
#endif
    MISC_SD_FREE(a_);
    MISC_SD_FREE(b_);
    MISC_SD_FREE(a_2);
    MISC_SD_FREE(b_2);


#endif
    /* --------------------------------------------------------- */






    /* --------------------------------------------------------------------- */
#if HAVE_MPI
FAILURE_1:
#endif
    MISC_SD_FREE(x);
    MISC_SD_FREE(y);
    MISC_SD_FREE(z);
    MISC_SD_FREE(ix);
    MISC_SD_FREE(iy);
    MISC_SD_FREE(iz);
    MISC_SD_FREE(ixy);
    MISC_SD_FREE(pnm0);
    MISC_SD_FREE(pnm1);
    MISC_SD_FREE(pnm2);
    MISC_SD_FREE(leg_cnm);
    MISC_SD_FREE(leg_snm);
    MISC_SD_FREE(ds);
#if DLAT > 0
    MISC_SD_FREE(tu);
    MISC_SD_FREE(dpnm0);
    MISC_SD_FREE(dpnm1);
    MISC_SD_FREE(dpnm2);
#endif
#if DLAT > 1
    MISC_SD_FREE(ddpnm0);
    MISC_SD_FREE(ddpnm1);
    MISC_SD_FREE(ddpnm2);
#endif
#if (DLAT > 0) || (DLON > 0)
    MISC_SD_FREE(u_rec);
#endif
#if (DLAT > 1) || (DLON > 1)
    MISC_SD_FREE(u2_rec);
#endif
#if KERNEL_GRAD > 0
    MISC_SD_FREE(leg_cnm_r);
    MISC_SD_FREE(leg_snm_r);
    MISC_SD_FREE(leg_cnm_p);
    MISC_SD_FREE(leg_snm_p);
#endif
#if KERNEL_GRAD > 1
    MISC_SD_FREE(leg_cnm_rr);
    MISC_SD_FREE(leg_snm_rr);
    MISC_SD_FREE(leg_cnm_rp);
    MISC_SD_FREE(leg_snm_rp);
    MISC_SD_FREE(leg_cnm_pp);
    MISC_SD_FREE(leg_snm_pp);
#endif
    /* --------------------------------------------------------------------- */






    return;
}
#endif

