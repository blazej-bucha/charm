/* This header file is not a part of API.
 *
 * Defines macros to compute sectorial, semi-sectorial and tesseral Legendre
 * functions.  We prefer macros over functions in this case to make sure the
 * code is always inlined in order to reduce the overhead.
 *
 * */


#ifndef __LEG_FUNC_XNUM_H__
#define __LEG_FUNC_XNUM_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


/* Macro to compute sectorial Legendre function. */
#define PNM_SECTORIAL_XNUM(x, ix, psmm1, ipsmm1, pnm0)                    \
                                                                          \
                                                                          \
         (x) = (psmm1);                                                   \
        (ix) = (ipsmm1);                                                  \
                                                                          \
                                                                          \
        if ((ix) == 0)                                                    \
            (pnm0) = (x);                                                 \
        else if ((ix) < -1)                                               \
            (pnm0) = ADDP(0.0);                                           \
        else if ((ix) < 0)                                                \
            (pnm0) = (x) * BIGI;                                          \
        else                                                              \
            (pnm0) = (x) * BIG;






/* Macro to compute semi-sectorial Legendre function. */
#define PNM_SEMISECTORIAL_XNUM(x, y, ix, iy, w, t, anm, pnm1)             \
                                                                          \
                                                                          \
         (y) =  (x);                                                      \
        (iy) = (ix);                                                      \
         (x) = ((anm) * (t)) * (y);                                       \
        (ix) = (iy);                                                      \
                                                                          \
                                                                          \
        (w) = FABS((x));                                                  \
        if ((w) >= BIGS)                                                  \
        {                                                                 \
             (x) *= BIGI;                                                 \
            (ix) += 1;                                                    \
        }                                                                 \
        else if ((w) < BIGSI)                                             \
        {                                                                 \
             (x) *= BIG;                                                  \
            (ix) -= 1;                                                    \
        }                                                                 \
                                                                          \
                                                                          \
        if ((ix) == 0)                                                    \
            (pnm1) = (x);                                                 \
        else if ((ix) < -1)                                               \
            (pnm1) = ADDP(0.0);                                           \
        else if ((ix) < 0)                                                \
            (pnm1) = (x) * BIGI;                                          \
        else                                                              \
            (pnm1) = (x) * BIG;






/* Macro to compute tesseral Legendre function.  Typical usage is inside the
 * innermost loop over harmonic degrees.
 *
 *      PNM_TESSERAL_XNUM(x, y, z, ix, iy, iz, ixy, w, t, anm[n], bnm[n], pnm2,
 *                        continue)
 *
 * The computed tesseral Legendre function will be stored in "pnm2".
 *
 * The last parameter called "UNDERFLOW" (see below) specifies what to do in
 * case the output of the macro, "pnm2", underflows.  In some cases, it may be
 * safe to use "continue", so that a few unnecessary multiplications by the
 * zero-valued "pnm2" are skipped.  These are usually present after the macro
 * is used in the main code.  In other cases, one can use "pnm2 = ADDP(0.0)",
 * where "pnm2" is the name of the variable that is supposed to store the
 * tesseral Legendre function in the main code.  The program will then evaluate
 * the entire macro and subsequently it will continue with commands found after
 * the macro.
 *
 * */
#define PNM_TESSERAL_XNUM(x, y, z, ix, iy, iz, ixy, w, t,                 \
                          anm, bnm, pnm2, UNDERFLOW)                      \
                                                                          \
                                                                          \
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
            UNDERFLOW;                                                    \
        else if ((iz) < 0)                                                \
            (pnm2) = (z) * BIGI;                                          \
        else                                                              \
            (pnm2) = (z) * BIG;


#ifdef __cplusplus
}
#endif


#endif

