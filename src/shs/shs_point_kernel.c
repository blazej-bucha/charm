/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "shs_check_symm_simd.h"
#include "../leg/leg_func_xnum.h"
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Returns summations over harmonic degrees, "a", "b", "a2" and "b2", for
 * synthesis with point values */
void CHARM(shs_point_kernel)(unsigned long nmax, unsigned long m,
                             const CHARM(shc) *shcs,
                             const REAL *anm, const REAL *bnm,
                             REAL_SIMD t, const REAL *ps, const int *ips,
                             REAL_SIMD ratio, REAL_SIMD ratio2,
                             REAL_SIMD ratiom, REAL_SIMD ratio2m,
                             REAL_SIMD symm_simd,
                             REAL_SIMD *a, REAL_SIMD *b,
                             REAL_SIMD *a2, REAL_SIMD *b2)
{
    REAL_SIMD w, x, y, z;
    RI_SIMD ixy, ix, iy, iz;
#ifdef SIMD
    RI_SIMD    zero_ri = SET_ZERO_RI;
    RI_SIMD    one_ri  = SET1_RI(1);
    RI_SIMD    mone_ri = SET1_RI(-1);
    REAL_SIMD  zero_r  = SET_ZERO_R;
    REAL_SIMD  BIG_r   = SET1_R(BIG);
    REAL_SIMD  BIGI_r  = SET1_R(BIGI);
    REAL_SIMD  BIGS_r  = SET1_R(BIGS);
    REAL_SIMD  BIGSI_r = SET1_R(BIGSI);
    REAL_SIMD  tmp1_r,  tmp2_r;
    MASK_SIMD  mask1, mask2;
    MASK2_SIMD mask3;
    ABS_R_INIT;
#endif


    REAL_SIMD pnm0, pnm1, pnm2;
    REAL_SIMD pnm_cnm, pnm_snm;
    REAL_SIMD ROOT3_r = SET1_R(ROOT3);


    REAL_SIMD ration  = ratiom;
    REAL_SIMD ratio2n = ratio2m;


    _Bool npm_even; /* True if "n + m" is even */
    unsigned long nmm; /* "n - m" */
    _Bool symm = CHARM(shs_check_symm_simd)(symm_simd);
    _Bool ds; /* Dynamical switching */


    /* Summations over harmonic degree "n" */
    if (m == 0)
    {
        /* Zonal harmonics */
        /* ----------------------------------------------------- */
        /* P00 */
        pnm0    = SET1_R(PREC(1.0));
        pnm_cnm = MUL_R(pnm0, SET1_R(shcs->c[0][0]));
        *a      = MUL_R(ration, pnm_cnm);
        *b      = SET_ZERO_R;
        ration  = MUL_R(ration, ratio);


        if (symm)
        {
            *a2     = MUL_R(ratio2n, pnm_cnm);
            *b2     = SET_ZERO_R;
            ratio2n = MUL_R(ratio2n, ratio2);
        }


        /* P10 */
        if (nmax >= 1)
        {
            pnm1    = MUL_R(ROOT3_r, t);
            pnm_cnm = MUL_R(pnm1, SET1_R(shcs->c[0][1]));
            *a      = ADD_R(*a, MUL_R(ration, pnm_cnm));
            ration  = MUL_R(ration, ratio);


            if (symm)
            {
                *a2     = SUB_R(*a2, MUL_R(ratio2n, pnm_cnm));
                ratio2n = MUL_R(ratio2n, ratio2);
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
                pnm2    = SUB_R(MUL_R(SET1_R(anm[n]), MUL_R(t, pnm1)),
                                MUL_R(SET1_R(bnm[n]), pnm0));
                pnm_cnm = MUL_R(pnm2, SET1_R(shcs->c[0][n]));
                *a      = ADD_R(*a, MUL_R(ration, pnm_cnm));
                ration  = MUL_R(ration, ratio);


                if (symm)
                {
                    if (npm_even)
                        *a2 = ADD_R(*a2, MUL_R(ratio2n, pnm_cnm));
                    else
                        *a2 = SUB_R(*a2, MUL_R(ratio2n, pnm_cnm));


                    ratio2n = MUL_R(ratio2n, ratio2);
                }


                pnm0 = pnm1;
                pnm1 = pnm2;
            }
        }
        /* ----------------------------------------------------- */
    }
    else /* Non-zonal harmonics */
    {

        /* Sectorial harmonics */
        /* ----------------------------------------------------- */
#ifdef SIMD
        PNM_SECTORIAL_XNUM_SIMD(x, ix, ps[(m - 1) * SIMD_SIZE],
                                ips[(m - 1) * SIMD_SIZE], pnm0,
                                BIG_r, zero_r, zero_ri, mone_ri, mask1, mask2,
                                SECTORIALS);
#else
        PNM_SECTORIAL_XNUM(x, ix, ps[m - 1], ips[m - 1], pnm0);
#endif


        pnm_cnm = MUL_R(pnm0, SET1_R(shcs->c[m][0]));
        pnm_snm = MUL_R(pnm0, SET1_R(shcs->s[m][0]));


        *a     = MUL_R(ration, pnm_cnm);
        *b     = MUL_R(ration, pnm_snm);
        ration = MUL_R(ration, ratio);


        if (symm)
        {
            *a2     = MUL_R(ratio2n, pnm_cnm);
            *b2     = MUL_R(ratio2n, pnm_snm);
            ratio2n = MUL_R(ratio2n, ratio2);
        }
        /* ----------------------------------------------------- */


        /* Tesseral harmonics */
        /* ----------------------------------------------------- */
        if (m < nmax)
        {
#ifdef SIMD
            PNM_SEMISECTORIAL_XNUM_SIMD(x, y, ix, iy, w, t, anm[m + 1],
                                        pnm1, mask1, mask2, mask3,
                                        zero_r, zero_ri, mone_ri, BIG_r,
                                        BIGS_r,  BIGI_r, SEMISECTORIALS);
#else
            PNM_SEMISECTORIAL_XNUM(x, y, ix, iy, w, t, anm[m + 1], pnm1);
#endif


            pnm_cnm = MUL_R(pnm1, SET1_R(shcs->c[m][1]));
            pnm_snm = MUL_R(pnm1, SET1_R(shcs->s[m][1]));


            *a     = ADD_R(*a, MUL_R(ration, pnm_cnm));
            *b     = ADD_R(*b, MUL_R(ration, pnm_snm));
            ration = MUL_R(ration, ratio);


            if (symm)
            {
                *a2     = SUB_R(*a2, MUL_R(ratio2n, pnm_cnm));
                *b2     = SUB_R(*b2, MUL_R(ratio2n, pnm_snm));
                ratio2n = MUL_R(ratio2n, ratio2);
            }


            /* Loop over degrees */
            /* ------------------------------------------------- */
            ds = 0;


            /* Is "n + m" even?  Since we start the loop with "n = m + 2", then
             * the parity of the first "m + 2 + m" is always even.  Then, it
             * changes with every loop iteration. */
            npm_even = 1;


            for (unsigned long n = (m + 2); n <= nmax;
                 n++, npm_even = !npm_even)
            {
                /* Compute tesseral Legendre function */
#ifdef SIMD
                PNM_TESSERAL_XNUM_SIMD(x, y, z, ix, iy, iz, ixy,
                                       w, t, anm[n], bnm[n], pnm2,
                                       tmp1_r, tmp2_r, mask1, mask2,
                                       mask3, zero_r, zero_ri, one_ri,
                                       BIG_r, BIGI_r, BIGS_r, BIGSI_r,
                                       TESSERALS1, TESSERALS2, ds);
#else
                PNM_TESSERAL_XNUM(x, y, z, ix, iy, iz,
                                  ixy, w, t, anm[n], bnm[n],
                                  pnm2, continue, ds);
#endif


                nmm = n - m;
                pnm_cnm = MUL_R(pnm2, SET1_R(shcs->c[m][nmm]));
                pnm_snm = MUL_R(pnm2, SET1_R(shcs->s[m][nmm]));


                *a     = ADD_R(*a, MUL_R(ration, pnm_cnm));
                *b     = ADD_R(*b, MUL_R(ration, pnm_snm));
                ration = MUL_R(ration, ratio);


                if (symm)
                {
                    if (npm_even)
                    {
                        *a2 = ADD_R(*a2, MUL_R(ratio2n, pnm_cnm));
                        *b2 = ADD_R(*b2, MUL_R(ratio2n, pnm_snm));
                    }
                    else
                    {
                        *a2 = SUB_R(*a2, MUL_R(ratio2n, pnm_cnm));
                        *b2 = SUB_R(*b2, MUL_R(ratio2n, pnm_snm));
                    }


                    ratio2n = MUL_R(ratio2n, ratio2);
                }


            } /* End of the loop over harmonic degrees */
            /* ------------------------------------------------- */


        } /* End of computation of tesseral harmonics */
        /* ----------------------------------------------------- */


    } /* End of the summation over harmonic degree "n" */
    /* --------------------------------------------------------- */






    return;
}

