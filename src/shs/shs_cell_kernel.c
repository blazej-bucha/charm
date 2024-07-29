/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "../leg/leg_func_xnum.h"
#include "../prec.h"
#include "../simd/simd.h"
#include "shs_check_symm_simd.h"
#include "shs_cell_kernel.h"
/* ------------------------------------------------------------------------- */






/* Returns summations over harmonic degrees, "a", "b", "a2" and "b2", for
 * synthesis with mean values.  The function also updates "imm0", "imm1" and
 * "imm2". */
void CHARM(shs_cell_kernel)(unsigned long nmax,
                            unsigned long m,
                            const CHARM(shc) *shcs,
                            const REAL *anm,
                            const REAL *bnm,
                            REAL_SIMD latmin,
                            REAL_SIMD latmax,
                            REAL_SIMD t1,
                            REAL_SIMD t2,
                            REAL_SIMD u1,
                            REAL_SIMD u2,
                            const REAL *ps1,
                            const REAL *ps2,
                            const int64_t *ips1,
                            const int64_t *ips2,
                            REAL_SIMD *imm0,
                            REAL_SIMD *imm1,
                            REAL_SIMD *imm2,
                            const REAL *en,
                            const REAL *fn,
                            const REAL *gm,
                            const REAL *hm,
                            const REAL *ri,
                            REAL_SIMD ratio,
                            REAL_SIMD ratio2,
                            REAL_SIMD ratiom,
                            REAL_SIMD ratio2m,
                            REAL_SIMD symm_simd,
                            REAL_SIMD *a,
                            REAL_SIMD *b,
                            REAL_SIMD *a2,
                            REAL_SIMD *b2)
{
    REAL_SIMD w;
    REAL_SIMD x1, x2, y1, y2, z1, z2;
    RI_SIMD ix1, ix2, iy1, iy2, iz1, iz2, ixy1, ixy2;
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
    NEG_R_INIT;
#endif


    REAL_SIMD pnm0_latmin, pnm1_latmin, pnm2_latmin;
    REAL_SIMD pnm0_latmax, pnm1_latmax, pnm2_latmax;
    REAL_SIMD in0, inm0, inm1, inm2;
    REAL_SIMD inm_cnm, inm_snm;


    REAL_SIMD ration  = ratiom;
    REAL_SIMD ratio2n = ratio2m;


    _Bool npm_even; /* True if "n + m" is even */
    unsigned long nmm; /* "n - m" */
    _Bool symm = CHARM(shs_check_symm_simd)(symm_simd);
    _Bool ds1, ds2; /* Dynamical switching */


    REAL_SIMD anms, bnms;


    iy1 = iy2 = iz1 = iz2 = ixy1 = ixy2 = SET_ZERO_RI;


    /* Computation of the lumped coefficients */
    if (m == 0)
    {


        /* Zonal Legendre functions and their integrals */
        /* ----------------------------------------------------------------- */
        /* P00 */
        /* ................................................................. */
        pnm0_latmin = SET1_R(PREC(1.0));
        pnm0_latmax = SET1_R(PREC(1.0));
        /* ................................................................. */


        /* P10 */
        /* ................................................................. */
        pnm1_latmin = t1;
        pnm1_latmax = t2;
        /* ................................................................. */


        /* I00 */
        /* ................................................................. */
        in0 = SUB_R(t2, t1);
        /* ................................................................. */


        /* Lumped coefficients */
        /* ................................................................. */
        inm_cnm = MUL_R(in0, SET1_R(shcs->c[0][0]));
        *a      = MUL_R(ration, inm_cnm);
        *b      = SET_ZERO_R;
        ration  = MUL_R(ration, ratio);


        if (symm)
        {
            *a2 = MUL_R(ratio2n, inm_cnm);
            *b2 = SET_ZERO_R;
            ratio2n = MUL_R(ratio2n, ratio2);
        }
        /* ................................................................. */



        if ((nmax + 1) >= 2)
        {
            /* is "n + m" even?  since we start the loop with "n = 1" and "m
             * = 0", then the parity of the first "n + m" is always odd.  then,
             * it changes with every loop iteration. */
            npm_even = 0;


            for (unsigned long n = 1; n <= nmax; n++, npm_even = !npm_even)
            {
                /* P20, P30, ..., Pnmax+1,0
                 * Since we use locally the "n" variable to compute Legendre
                 * polynomials of degree "n + 1", the "n + 1"th elements have
                 * to be taken from the vectors "en" and "fn" below. Once
                 * again, note that when "m == 0", then "pnm0", "pnm1" and
                 * "pnm2" for latmin and latmax represent un-normalized
                 * Legendre polynomials. These are needed to get the integrals
                 * of fully-normalized Legendre polynomials */
                /* ......................................................... */
                pnm2_latmin = SUB_R(MUL_R(SET1_R(en[n + 1]),
                                          MUL_R(t1, pnm1_latmin)),
                                    MUL_R(SET1_R(fn[n + 1]), pnm0_latmin));
                pnm2_latmax = SUB_R(MUL_R(SET1_R(en[n + 1]),
                                          MUL_R(t2, pnm1_latmax)),
                                    MUL_R(SET1_R(fn[n + 1]), pnm0_latmax));
                /* ......................................................... */


                /* I10, I20, ..., Inmax,0
                 * Computed from Pn+1,0 and Pn-1,0 */
                /* ......................................................... */
                in0 = MUL_R(SET1_R(ri[2 * n + 1]),
                            SUB_R(SUB_R(pnm2_latmax, pnm0_latmax),
                                  SUB_R(pnm2_latmin, pnm0_latmin)));
                /* ......................................................... */


                /* Lumped coefficients */
                /* ......................................................... */
                inm_cnm = MUL_R(in0, SET1_R(shcs->c[0][n]));
                *a      = ADD_R(*a, MUL_R(ration, inm_cnm));
                ration  = MUL_R(ration, ratio);

                if (symm)
                {
                    if (npm_even)
                        *a2 = ADD_R(*a2, MUL_R(ratio2n, inm_cnm));
                    else
                        *a2 = SUB_R(*a2, MUL_R(ratio2n, inm_cnm));


                    ratio2n = MUL_R(ratio2n, ratio2);
                }
                /* ......................................................... */


                pnm0_latmin = pnm1_latmin;
                pnm1_latmin = pnm2_latmin;


                pnm0_latmax = pnm1_latmax;
                pnm1_latmax = pnm2_latmax;

            }
        }
        /* ----------------------------------------------------------------- */


    }
    else /* Non-zonal harmonics */
    {


        /* Sectorial Legendre functions and their integrals */
        /* ----------------------------------------------------------------- */

        /* Pmm for "latmin" */
#ifdef SIMD
        PNM_SECTORIAL_XNUM_SIMD(x1, ix1, ps1[(m - 1) * SIMD_SIZE],
                                ips1[(m - 1) * SIMD_SIZE], pnm0_latmin,
                                BIG_r, zero_r, zero_ri, mone_ri, mask1, mask2,
                                SECTORIALS1);
#else
        PNM_SECTORIAL_XNUM(x1, ix1, ps1[m - 1], ips1[m - 1], pnm0_latmin);
#endif


        /* Pmm for "latmax" */
#ifdef SIMD
        PNM_SECTORIAL_XNUM_SIMD(x2, ix2, ps2[(m - 1) * SIMD_SIZE],
                                ips2[(m - 1) * SIMD_SIZE], pnm0_latmax,
                                BIG_r, zero_r, zero_ri, mone_ri, mask1, mask2,
                                SECTORIALS2);
#else
        PNM_SECTORIAL_XNUM(x2, ix2, ps2[m - 1], ips2[m - 1], pnm0_latmax);
#endif


        /* Imm */
        /* ................................................................. */
        if (m == 1) /* I11 */
        {
            *imm0 = MUL_R(SET1_R(ROOT3 / PREC(2.0)),
                          SUB_R(SUB_R(MUL_R(t2, u2),
                                      SUB_R(SET1_R(PI_2), latmax)),
                                SUB_R(MUL_R(t1, u1),
                                      SUB_R(SET1_R(PI_2), latmin))));
            inm0 = *imm0;
        }
        else if (m == 2) /* I22 */
        {
            *imm1 = MUL_R(SET1_R(SQRT(PREC(15.0)) / PREC(6.0)),
                          SUB_R(MUL_R(t2, SUB_R(SET1_R(PREC(3.0)),
                                                MUL_R(t2, t2))),
                                MUL_R(t1, SUB_R(SET1_R(PREC(3.0)),
                                                MUL_R(t1, t1)))));
            inm0 = *imm1;
        }
        else /* I33, I44, ..., Inmax,nmax */
        {
            *imm2 = ADD_R(MUL_R(SET1_R(gm[m]), *imm0),
                          MUL_R(SET1_R(PREC(1.0) / (REAL)(m + 1)),
                                SUB_R(MUL_R(t2, pnm0_latmax),
                                      MUL_R(t1, pnm0_latmin))));
            inm0 = *imm2;

            *imm0 = *imm1;
            *imm1 = *imm2;
        }
        /* ................................................................. */


        /* Lumped coefficients */
        /* ................................................................. */
        inm_cnm = MUL_R(inm0, SET1_R(shcs->c[m][0]));
        inm_snm = MUL_R(inm0, SET1_R(shcs->s[m][0]));


        *a     = MUL_R(ration, inm_cnm);
        *b     = MUL_R(ration, inm_snm);
        ration = MUL_R(ration, ratio);


        if (symm)
        {
            *a2     = MUL_R(ratio2n, inm_cnm);
            *b2     = MUL_R(ratio2n, inm_snm);
            ratio2n = MUL_R(ratio2n, ratio2);
        }
        /* ................................................................. */
        /* ----------------------------------------------------------------- */






        /* Tesseral harmonics */
        /* ----------------------------------------------------------------- */
        if (m < nmax)
        {
            anms = SET1_R(anm[m + 1]);
            bnms = SET1_R(bnm[m + 1]);


            /* Pm+1,m for "latmin" */
#ifdef SIMD
            PNM_SEMISECTORIAL_XNUM_SIMD(x1, y1, ix1, iy1, w, t1, anms,
                                        pnm1_latmin, mask1, mask2, mask3,
                                        zero_r, zero_ri, mone_ri, BIG_r,
                                        BIGS_r,  BIGI_r, SEMISECTORIALS1);
#else
            PNM_SEMISECTORIAL_XNUM(x1, y1, ix1, iy1, w, t1, anms,
                                   pnm1_latmin);
#endif


            /* Pm+1,m for "latmax" */
#ifdef SIMD
            PNM_SEMISECTORIAL_XNUM_SIMD(x2, y2, ix2, iy2, w, t2, anms,
                                        pnm1_latmax, mask1, mask2, mask3,
                                        zero_r, zero_ri, mone_ri, BIG_r,
                                        BIGS_r,  BIGI_r, SEMISECTORIALS2);
#else
            PNM_SEMISECTORIAL_XNUM(x2, y2, ix2, iy2, w, t2, anms,
                                   pnm1_latmax);
#endif


            /* Im+1,m */
            /* ............................................................. */
            /* This is not a typo, "pnm0" are indeed required here */
            inm1 = NEG_R(MUL_R(SET1_R(anm[m + 1] / (REAL)(m + 2)),
                               SUB_R(MUL_R(MUL_R(u2, u2), pnm0_latmax),
                                     MUL_R(MUL_R(u1, u1),
                                           pnm0_latmin))));
            /* ............................................................. */


            /* Lumped coefficients */
            /* ............................................................. */
            inm_cnm = MUL_R(inm1, SET1_R(shcs->c[m][1]));
            inm_snm = MUL_R(inm1, SET1_R(shcs->s[m][1]));


            *a     = ADD_R(*a, MUL_R(ration, inm_cnm));
            *b     = ADD_R(*b, MUL_R(ration, inm_snm));
            ration = MUL_R(ration, ratio);


            if (symm)
            {
                *a2     = SUB_R(*a2, MUL_R(ratio2n, inm_cnm));
                *b2     = SUB_R(*b2, MUL_R(ratio2n, inm_snm));
                ratio2n = MUL_R(ratio2n, ratio2);
            }
            /* ............................................................. */


            /* Pm+2,m, Pm+3,m, ..., Pnmax,m and their integrals */
            /* ............................................................. */
            ds1 = ds2 = 0;


            /* Is "n + m" even?  Since we start the loop with "n = m + 2", then
             * the parity of the first "m + 2 + m" is always even.  Then, it
             * changes with every loop iteration. */
            npm_even = 1;


            for (unsigned long n = (m + 2); n <= nmax;
                 n++, npm_even = !npm_even)
            {
                anms = SET1_R(anm[n]);
                bnms = SET1_R(bnm[n]);


                /* Pm+2,m, Pm+3,m, ..., Pnmax,m for "latmin" */
#ifdef SIMD
                PNM_TESSERAL_XNUM_SIMD(x1, y1, z1, ix1, iy1, iz1, ixy1,
                                       w, t1, anms, bnms, pnm2_latmin,
                                       tmp1_r, tmp2_r, mask1, mask2,
                                       mask3, zero_r, zero_ri, one_ri,
                                       BIG_r, BIGI_r, BIGS_r, BIGSI_r,
                                       TESSERALS1, TESSERALS2, ds1);
#else
                PNM_TESSERAL_XNUM(x1, y1, z1, ix1, iy1, iz1,
                                  ixy1, w, t1, anms, bnms,
                                  pnm2_latmin, ds1);
#endif


                /* Pm+2,m, Pm+3,m, ..., Pnmax,m for "latmax" */
#ifdef SIMD
                PNM_TESSERAL_XNUM_SIMD(x2, y2, z2, ix2, iy2, iz2, ixy2,
                                       w, t2, anms, bnms, pnm2_latmax,
                                       tmp1_r, tmp2_r, mask1, mask2,
                                       mask3, zero_r, zero_ri, one_ri,
                                       BIG_r, BIGI_r, BIGS_r, BIGSI_r,
                                       TESSERALS3, TESSERALS4, ds2);
#else
                PNM_TESSERAL_XNUM(x2, y2, z2, ix2, iy2, iz2,
                                  ixy2, w, t2, anms, bnms,
                                  pnm2_latmax, ds2);
#endif


                /* Im+2,m, Im+3,m, ..., Inmax,m */
                /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                /* This is not a typo, "pnm1" are indeed required here */
                inm2 = SUB_R(MUL_R(MUL_R(SET1_R(hm[n]),
                                         bnms), inm0),
                             MUL_R(SET1_R(anm[n] / (REAL)(n + 1)),
                                   SUB_R(MUL_R(MUL_R(u2, u2),
                                               pnm1_latmax),
                                         MUL_R(MUL_R(u1, u1),
                                               pnm1_latmin))));

                inm0 = inm1;
                inm1 = inm2;
                /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


                /* Lumped coefficients */
                /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                nmm = n - m;
                inm_cnm = MUL_R(inm2, SET1_R(shcs->c[m][nmm]));
                inm_snm = MUL_R(inm2, SET1_R(shcs->s[m][nmm]));


                *a     = ADD_R(*a, MUL_R(ration, inm_cnm));
                *b     = ADD_R(*b, MUL_R(ration, inm_snm));
                ration = MUL_R(ration, ratio);


                if (symm)
                {
                    if (npm_even)
                    {
                        *a2 = ADD_R(*a2, MUL_R(ratio2n, inm_cnm));
                        *b2 = ADD_R(*b2, MUL_R(ratio2n, inm_snm));
                    }
                    else
                    {
                        *a2 = SUB_R(*a2, MUL_R(ratio2n, inm_cnm));
                        *b2 = SUB_R(*b2, MUL_R(ratio2n, inm_snm));
                    }


                    ratio2n = MUL_R(ratio2n, ratio2);
                }
                /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


                pnm1_latmin = pnm2_latmin;
                pnm1_latmax = pnm2_latmax;


            } /* End of the loop over harmonic degrees */
            /* ............................................................. */


        } /* End of computation of tesseral harmonics */
        /* ----------------------------------------------------------------- */


    } /* End of computation of the lumped coefficients */
    /* --------------------------------------------------------------------- */






    return;
}

