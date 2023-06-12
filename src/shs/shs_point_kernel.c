/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shs_check_symm_simd.h"
#include "../leg/leg_func_xnum.h"
#include "../leg/leg_func_use_xnum.h"
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
/* The "LOOP_ITER" macro can be used for any values of "r / R".
 *
 * The "LOOP_ITER_R1" macro can *only* be used if all "(r / R) == 1.0".
 *
 * "LOOP_ITER_R1" is faster than "LOOP_ITER", so the "is_ratio_one" value is
 * used to check whether "LOOP_ITER_R1" can safely be used or not.
 *
 * Except for "LOOP_ITER_R1", the entire kernel assumes that "(r / R) != 1", so
 * that the slower code variant is used elsewhere.  This include code to
 * compute Legendre polynomials, sectorial Legendre functions, and tesseral
 * Legendre functions (the last one with X-numbers).  This is because the
 * computing time of these parts of the code are rather negligible and using
 * a separate code for "(r / R) == 1.0" would make the code uselessly
 * lengthy. */


#define LOOP_ITER(n, PM_R)                                                    \
    anms = SET1_R(anm[(n)]);                                                  \
    bnms = SET1_R(bnm[(n)]);                                                  \
    cnm  = SET1_R(shcs->c[m][idx]);                                           \
    snm  = SET1_R(shcs->s[m][idx]);                                           \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        PNM_TESSERAL(x[l], y[l], pnm2[l], t[l], anms, bnms);                  \
    }                                                                         \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        pnm_cnm[l] = MUL_R(pnm2[l], cnm);                                     \
    }                                                                         \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        pnm_snm[l] = MUL_R(pnm2[l], snm);                                     \
    }                                                                         \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        a[l] = ADD_R(a[l], MUL_R(ration[l], pnm_cnm[l]));                     \
    }                                                                         \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        b[l] = ADD_R(b[l], MUL_R(ration[l], pnm_snm[l]));                     \
    }                                                                         \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        ration[l] = MUL_R(ration[l], ratio[l]);                               \
    }                                                                         \
                                                                              \
                                                                              \
    if (symm)                                                                 \
    {                                                                         \
        for (l = 0; l < SIMD_BLOCK; l++)                                      \
        {                                                                     \
            a2[l] = PM_R(a2[l], MUL_R(ratio2n[l], pnm_cnm[l]));               \
        }                                                                     \
        for (l = 0; l < SIMD_BLOCK; l++)                                      \
        {                                                                     \
            b2[l] = PM_R(b2[l], MUL_R(ratio2n[l], pnm_snm[l]));               \
        }                                                                     \
                                                                              \
                                                                              \
        for (l = 0; l < SIMD_BLOCK; l++)                                      \
        {                                                                     \
            ratio2n[l] = MUL_R(ratio2n[l], ratio2[l]);                        \
        }                                                                     \
    }                                                                         \
                                                                              \
                                                                              \
    idx++;




#define LOOP_ITER_R1(n, PM_R)                                                 \
    anms = SET1_R(anm[(n)]);                                                  \
    bnms = SET1_R(bnm[(n)]);                                                  \
    cnm  = SET1_R(shcs->c[m][idx]);                                           \
    snm  = SET1_R(shcs->s[m][idx]);                                           \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        PNM_TESSERAL(x[l], y[l], pnm2[l], t[l], anms, bnms);                  \
    }                                                                         \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        pnm_cnm[l] = MUL_R(pnm2[l], cnm);                                     \
    }                                                                         \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        pnm_snm[l] = MUL_R(pnm2[l], snm);                                     \
    }                                                                         \
                                                                              \
                                                                              \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        a[l] = ADD_R(a[l], pnm_cnm[l]);                                       \
    }                                                                         \
    for (l = 0; l < SIMD_BLOCK; l++)                                          \
    {                                                                         \
        b[l] = ADD_R(b[l], pnm_snm[l]);                                       \
    }                                                                         \
                                                                              \
                                                                              \
    if (symm)                                                                 \
    {                                                                         \
        for (l = 0; l < SIMD_BLOCK; l++)                                      \
        {                                                                     \
            a2[l] = PM_R(a2[l], pnm_cnm[l]);                                  \
        }                                                                     \
        for (l = 0; l < SIMD_BLOCK; l++)                                      \
        {                                                                     \
            b2[l] = PM_R(b2[l], pnm_snm[l]);                                  \
        }                                                                     \
    }                                                                         \
                                                                              \
                                                                              \
    idx++;
/* ------------------------------------------------------------------------- */






/* Returns summations over harmonic degrees, "a", "b", "a2" and "b2", for
 * synthesis with point values */
void CHARM(shs_point_kernel)(unsigned long nmax, unsigned long m,
                             const CHARM(shc) *shcs, _Bool is_ratio_one,
                             const REAL *anm, const REAL *bnm,
                             REAL_SIMD *t, const REAL *ps, const int *ips,
                             REAL_SIMD *ratio, REAL_SIMD *ratio2,
                             REAL_SIMD *ratiom, REAL_SIMD *ratio2m,
                             REAL_SIMD *symm_simd,
                             REAL_SIMD *a, REAL_SIMD *b,
                             REAL_SIMD *a2, REAL_SIMD *b2)
{
    REAL_SIMD w, x[SIMD_BLOCK], y[SIMD_BLOCK], z[SIMD_BLOCK];
    RI_SIMD ixy[SIMD_BLOCK], ix[SIMD_BLOCK], iy[SIMD_BLOCK], iz[SIMD_BLOCK];
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
    ABS_R_INIT;
#endif


    REAL_SIMD pnm0[SIMD_BLOCK], pnm1[SIMD_BLOCK], pnm2[SIMD_BLOCK];
    REAL_SIMD pnm_cnm[SIMD_BLOCK], pnm_snm[SIMD_BLOCK];
    const REAL_SIMD ROOT3_r = SET1_R(ROOT3);


    size_t l;
    REAL_SIMD ration[SIMD_BLOCK];
    REAL_SIMD ratio2n[SIMD_BLOCK];


    for (l = 0; l < SIMD_BLOCK; l++)
    {
        ration[l]  = ratiom[l];
        ratio2n[l] = ratio2m[l];
    }


    REAL_SIMD cnm, snm;


    _Bool npm_even; /* True if "n + m" is even */
    _Bool symm = 0;
    for (l = 0; l < SIMD_BLOCK; l++)
        symm = symm || CHARM(shs_check_symm_simd)(symm_simd[l]);
    _Bool ds[SIMD_BLOCK]; /* Dynamical switching */


    REAL_SIMD anms, bnms;


    /* Summations over harmonic degree "n" */
    if (m == 0)
    {
        /* Zonal harmonics */
        /* ----------------------------------------------------- */
        /* P00 */
        cnm = SET1_R(shcs->c[0][0]);
        for (l = 0; l < SIMD_BLOCK; l++)
        {
            pnm0[l]    = SET1_R(PREC(1.0));
            pnm_cnm[l] = MUL_R(pnm0[l], cnm);
            a[l]       = MUL_R(ration[l], pnm_cnm[l]);
            b[l]       = SET_ZERO_R;
            ration[l]  = MUL_R(ration[l], ratio[l]);


            if (symm)
            {
                a2[l]      = MUL_R(ratio2n[l], pnm_cnm[l]);
                b2[l]      = SET_ZERO_R;
                ratio2n[l] = MUL_R(ratio2n[l], ratio2[l]);
            }
        }


        /* P10 */
        if (nmax >= 1)
        {
            cnm = SET1_R(shcs->c[0][1]);
            for (l = 0; l < SIMD_BLOCK; l++)
            {
                pnm1[l]    = MUL_R(ROOT3_r, t[l]);
                pnm_cnm[l] = MUL_R(pnm1[l], cnm);
                a[l]       = ADD_R(a[l], MUL_R(ration[l], pnm_cnm[l]));
                ration[l]  = MUL_R(ration[l], ratio[l]);


                if (symm)
                {
                    a2[l]      = SUB_R(a2[l], MUL_R(ratio2n[l], pnm_cnm[l]));
                    ratio2n[l] = MUL_R(ratio2n[l], ratio2[l]);
                }
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
                cnm  = SET1_R(shcs->c[0][n]);


                for (l = 0; l < SIMD_BLOCK; l++)
                {
                    pnm2[l]    = SUB_R(MUL_R(anms,
                                             MUL_R(t[l], pnm1[l])),
                                       MUL_R(bnms, pnm0[l]));
                    pnm_cnm[l] = MUL_R(pnm2[l], cnm);
                    a[l]       = ADD_R(a[l], MUL_R(ration[l], pnm_cnm[l]));
                    ration[l]  = MUL_R(ration[l], ratio[l]);


                    pnm0[l] = pnm1[l];
                    pnm1[l] = pnm2[l];


                    if (symm)
                    {
                        if (npm_even)
                            a2[l] = ADD_R(a2[l],
                                          MUL_R(ratio2n[l], pnm_cnm[l]));
                        else
                            a2[l] = SUB_R(a2[l],
                                          MUL_R(ratio2n[l], pnm_cnm[l]));


                        ratio2n[l] = MUL_R(ratio2n[l], ratio2[l]);
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
        cnm = SET1_R(shcs->c[m][0]);
        snm = SET1_R(shcs->s[m][0]);


        for (l = 0; l < SIMD_BLOCK; l++)
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


            pnm_cnm[l] = MUL_R(pnm0[l], cnm);
            pnm_snm[l] = MUL_R(pnm0[l], snm);


            a[l]      = MUL_R(ration[l], pnm_cnm[l]);
            b[l]      = MUL_R(ration[l], pnm_snm[l]);
            ration[l] = MUL_R(ration[l], ratio[l]);


            if (symm)
            {
                a2[l]      = MUL_R(ratio2n[l], pnm_cnm[l]);
                b2[l]      = MUL_R(ratio2n[l], pnm_snm[l]);
                ratio2n[l] = MUL_R(ratio2n[l], ratio2[l]);
            }
        }
        /* ----------------------------------------------------- */


        /* Tesseral harmonics */
        /* ----------------------------------------------------- */
        if (m < nmax)
        {
            anms = SET1_R(anm[m + 1]);
            bnms = SET1_R(bnm[m + 1]);
            cnm  = SET1_R(shcs->c[m][1]);
            snm  = SET1_R(shcs->s[m][1]);


            for (l = 0; l < SIMD_BLOCK; l++)
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


                pnm_cnm[l] = MUL_R(pnm1[l], cnm);
                pnm_snm[l] = MUL_R(pnm1[l], snm);


                a[l]      = ADD_R(a[l], MUL_R(ration[l], pnm_cnm[l]));
                b[l]      = ADD_R(b[l], MUL_R(ration[l], pnm_snm[l]));
                ration[l] = MUL_R(ration[l], ratio[l]);


                if (symm)
                {
                    a2[l]      = SUB_R(a2[l], MUL_R(ratio2n[l], pnm_cnm[l]));
                    b2[l]      = SUB_R(b2[l], MUL_R(ratio2n[l], pnm_snm[l]));
                    ratio2n[l] = MUL_R(ratio2n[l], ratio2[l]);
                }
            }


            /* Loop over degrees */
            /* ------------------------------------------------- */
            for (l = 0; l < SIMD_BLOCK; l++)
                ds[l] = 0;


            /* Is "n + m" even?  Since we start the loop with "n = m + 2", then
             * the parity of the first "m + 2 + m" is always even.  Then, it
             * changes with every loop iteration. */
            npm_even = 1;


            unsigned long n;
            unsigned long idx = 2;
            for (n = (m + 2);
                 CHARM(leg_func_use_xnum(ds, SIMD_BLOCK)) && n <= nmax;
                 n++, npm_even = !npm_even, idx++)
            {
                anms = SET1_R(anm[n]);
                bnms = SET1_R(bnm[n]);
                cnm  = SET1_R(shcs->c[m][idx]);
                snm  = SET1_R(shcs->s[m][idx]);


                /* Compute tesseral Legendre function */
                for (l = 0; l < SIMD_BLOCK; l++)
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


                for (l = 0; l < SIMD_BLOCK; l++)
                    pnm_cnm[l] = MUL_R(pnm2[l], cnm);
                for (l = 0; l < SIMD_BLOCK; l++)
                    pnm_snm[l] = MUL_R(pnm2[l], snm);


                for (l = 0; l < SIMD_BLOCK; l++)
                    a[l] = ADD_R(a[l], MUL_R(ration[l], pnm_cnm[l]));
                for (l = 0; l < SIMD_BLOCK; l++)
                    b[l] = ADD_R(b[l], MUL_R(ration[l], pnm_snm[l]));


                for (l = 0; l < SIMD_BLOCK; l++)
                    ration[l] = MUL_R(ration[l], ratio[l]);


                if (symm)
                {
                    if (npm_even)
                    {
                        for (l = 0; l < SIMD_BLOCK; l++)
                            a2[l] = ADD_R(a2[l],
                                          MUL_R(ratio2n[l], pnm_cnm[l]));
                        for (l = 0; l < SIMD_BLOCK; l++)
                            b2[l] = ADD_R(b2[l],
                                          MUL_R(ratio2n[l], pnm_snm[l]));
                    }
                    else
                    {
                        for (l = 0; l < SIMD_BLOCK; l++)
                            a2[l] = SUB_R(a2[l],
                                          MUL_R(ratio2n[l], pnm_cnm[l]));
                        for (l = 0; l < SIMD_BLOCK; l++)
                            b2[l] = SUB_R(b2[l],
                                          MUL_R(ratio2n[l], pnm_snm[l]));
                    }


                    for (l = 0; l < SIMD_BLOCK; l++)
                        ratio2n[l] = MUL_R(ratio2n[l], ratio2[l]);
                }
            }


            if (n > nmax)
                return;


            /* From now on, "F"-numbers can be used instead of the "X"-numbers
             * to gain some speed */


            /* We want to unroll the loop that follows.  To do that, we need to
             * make sure that the loop starts with an even value of "n + m".
             * So if "n + m" is odd, we need to hard code one iteration. */
            if (!npm_even)
            {
                LOOP_ITER(n, SUB_R);
                n++;
            }


            if (is_ratio_one)
            {
                for (; (n + 1) <= nmax; n += 2)
                {
                    LOOP_ITER_R1(n,     ADD_R);
                    LOOP_ITER_R1(n + 1, SUB_R);
                }
            }
            else
            {
                for (; (n + 1) <= nmax; n += 2)
                {
                    LOOP_ITER(n,     ADD_R);
                    LOOP_ITER(n + 1, SUB_R);
                }
            }


            if (n > nmax)
                return;


            LOOP_ITER(n, ADD_R);
            /* ------------------------------------------------- */


        } /* End of computation of tesseral harmonics */
        /* ----------------------------------------------------- */


    } /* End of the summation over harmonic degree "n" */
    /* --------------------------------------------------------- */






    return;
}

