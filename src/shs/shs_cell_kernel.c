/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../leg/leg_func_xnum.h"
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Returns summations over harmonic degrees, "a", "b", "a2" and "b2", for
 * synthesis with mean values.  The function also updates "imm0", "imm1" and
 * "imm2". */
void CHARM(shs_cell_kernel)(unsigned long nmax, unsigned long m,
                            const CHARM(shc) *shcs,
                            const REAL *anm, const REAL *bnm,
                            REAL latmini, REAL latmaxi,
                            REAL t1, REAL t2,
                            REAL u1, REAL u2,
                            const REAL *ps1, const REAL *ps2,
                            const int *ips1, const int *ips2,
                            REAL *imm0, REAL *imm1, REAL *imm2,
                            const REAL *en, const REAL *fn,
                            const REAL *gm, const REAL *hm,
                            const REAL *ri, const REAL *rpows,
                            const REAL *rpows2,
                            _Bool symmi,
                            REAL *a,  REAL *b,
                            REAL *a2, REAL *b2)
{
    REAL w;
    REAL x1, x2, y1, y2, z1, z2;
    int ix1, ix2, iy1, iy2, iz1, iz2, ixy1, ixy2;


    REAL pnm0_latmini, pnm1_latmini, pnm2_latmini;
    REAL pnm0_latmaxi, pnm1_latmaxi, pnm2_latmaxi;
    REAL in0, inm0, inm1, inm2;
    REAL inm_cnm, inm_snm;


    _Bool npm_even; /* True if "n + m" is even */


    iy1 = iy2 = iz1 = iz2 = ixy1 = ixy2 = 0;


    /* Computation of the lumped coefficients */
    if (m == 0)
    {


        /* Zonal Legendre functions and their integrals */
        /* ----------------------------------------------------------------- */
        /* P00 */
        /* ................................................................. */
        pnm0_latmini = ADDP(1.0);
        pnm0_latmaxi = ADDP(1.0);
        /* ................................................................. */


        /* P10 */
        /* ................................................................. */
        pnm1_latmini = t1;
        pnm1_latmaxi = t2;
        /* ................................................................. */


        /* I00 */
        /* ................................................................. */
        in0 = t2 - t1;
        /* ................................................................. */


        /* Lumped coefficients */
        /* ................................................................. */
        inm_cnm = in0 * shcs->c[0][0];
        *a = rpows[0] * inm_cnm;
        *b = ADDP(0.0);


        if (symmi)
        {
            *a2 = rpows2[0] * inm_cnm;
            *b2 = ADDP(0.0);
        }
        /* ................................................................. */



        if ((nmax + 1) >= 2)
        {
            for (unsigned long n = 1; n <= nmax; n++)
            {
                /* P20, P30, ..., Pnmax+1,0
                 * Since we use locally the "n" variable to compute Legendre
                 * polynomials of degree "n + 1", the "n + 1"th elements have
                 * to be taken from the vectors "en" and "fn" below. Once
                 * again, note that when "m == 0", then "pnm0", "pnm1" and
                 * "pnm2" for latmini and latmaxi represent un-normalized
                 * Legendre polynomials. These are needed to get the integrals
                 * of fully-normalized Legendre polynomials */
                /* ......................................................... */
                pnm2_latmini = en[n + 1] * t1 * pnm1_latmini
                             - fn[n + 1] * pnm0_latmini;
                pnm2_latmaxi = en[n + 1] * t2 * pnm1_latmaxi
                             - fn[n + 1] * pnm0_latmaxi;
                /* ......................................................... */


                /* I10, I20, ..., Inmax,0
                 * Computed from Pn+1,0 and Pn-1,0 */
                /* ......................................................... */
                in0 = ri[2 * n + 1] * (pnm2_latmaxi - pnm0_latmaxi -
                                       pnm2_latmini + pnm0_latmini);
                /* ......................................................... */


                /* Lumped coefficients */
                /* ......................................................... */
                inm_cnm = in0 * shcs->c[0][n];
                *a     += rpows[n] * inm_cnm;

                if (symmi)
                {
                    if ((n % 2) == 0)
                        *a2 += rpows2[n] * inm_cnm;
                    else
                        *a2 -= rpows2[n] * inm_cnm;
                }
                /* ......................................................... */


                pnm0_latmini = pnm1_latmini;
                pnm1_latmini = pnm2_latmini;


                pnm0_latmaxi = pnm1_latmaxi;
                pnm1_latmaxi = pnm2_latmaxi;

            }
        }
        /* ----------------------------------------------------------------- */


    }
    else /* Non-zonal harmonics */
    {


        /* Sectorial Legendre functions and their integrals */
        /* ----------------------------------------------------------------- */

        /* Pmm for latmini */
        PNM_SECTORIAL_XNUM(x1, ix1, ps1[m - 1], ips1[m - 1], pnm0_latmini);


        /* Pmm for latmaxi */
        PNM_SECTORIAL_XNUM(x2, ix2, ps2[m - 1], ips2[m - 1], pnm0_latmaxi);


        /* Imm */
        /* ................................................................. */
        if (m == 1) /* I11 */
        {
            *imm0 = SQRT(ADDP(3.0)) / ADDP(2.0) *
                         ((t2 * u2 - (PI_2 - latmaxi)) -
                          (t1 * u1 - (PI_2 - latmini)));
            inm0 = *imm0;
        }
        else if (m == 2) /* I22 */
        {
            *imm1 = SQRT(ADDP(15.0)) / ADDP(6.0) *
                         (t2 * (ADDP(3.0) - t2 * t2) -
                          t1 * (ADDP(3.0) - t1 * t1));
            inm0 = *imm1;
        }
        else /* I33, I44, ..., Inmax,nmax */
        {
            *imm2 = gm[m] * (*imm0) + ADDP(1.0) / (REAL)(m + 1) *
                    (t2 * pnm0_latmaxi - t1 * pnm0_latmini);
            inm0 = *imm2;

            *imm0 = *imm1;
            *imm1 = *imm2;
        }
        /* ................................................................. */


        /* Lumped coefficients */
        /* ................................................................. */
        inm_cnm = inm0 * shcs->c[m][0];
        inm_snm = inm0 * shcs->s[m][0];


        *a = rpows[m] * inm_cnm;
        *b = rpows[m] * inm_snm;


        if (symmi)
        {
            *a2 = rpows2[m] * inm_cnm;
            *b2 = rpows2[m] * inm_snm;
        }
        /* ................................................................. */
        /* ----------------------------------------------------------------- */






        /* Tesseral harmonics */
        /* ----------------------------------------------------------------- */
        if (m < nmax)
        {

            /* Pm+1,m for latmini */
            PNM_SEMISECTORIAL_XNUM(x1, y1, ix1, iy1, w, t1, anm[m + 1],
                                   pnm1_latmini);


            /* Pm+1,m for latmaxi */
            PNM_SEMISECTORIAL_XNUM(x2, y2, ix2, iy2, w, t2, anm[m + 1],
                                   pnm1_latmaxi);


            /* Im+1,m */
            /* ............................................................. */
            /* This is not a typo, "pnm0" are indeed required here */
            inm1 = -(anm[m + 1] / (REAL)(m + 2)) *
                    ((u2 * u2) * pnm0_latmaxi - (u1 * u1) * pnm0_latmini);
            /* ............................................................. */


            /* Lumped coefficients */
            /* ............................................................. */
            inm_cnm = inm1 * shcs->c[m][1];
            inm_snm = inm1 * shcs->s[m][1];


            *a += rpows[m + 1] * inm_cnm;
            *b += rpows[m + 1] * inm_snm;


            if (symmi)
            {
                *a2 -= rpows2[m + 1] * inm_cnm;
                *b2 -= rpows2[m + 1] * inm_snm;
            }
            /* ............................................................. */


            /* Pm+2,m, Pm+3,m, ..., Pnmax,m and their integrals */
            /* ............................................................. */
            /* Is "n + m" even?  Since we start the loop with "n = m + 2", then
             * the parity of the first "m + 2 + m" is always even.  Then, it
             * changes with every loop iteration. */
            npm_even = 1;
            for (unsigned long n = (m + 2); n <= nmax;
                 n++, npm_even = !npm_even)
            {

                /* Pm+2,m, Pm+3,m, ..., Pnmax,m for latmini */
                PNM_TESSERAL_XNUM(x1, y1, z1, ix1, iy1, iz1,
                                  ixy1, w, t1, anm[n], bnm[n],
                                  pnm2_latmini, pnm2_latmini = ADDP(0.0));


                /* Pm+2,m, Pm+3,m, ..., Pnmax,m for latmaxi */
                PNM_TESSERAL_XNUM(x2, y2, z2, ix2, iy2, iz2,
                                  ixy2, w, t2, anm[n], bnm[n],
                                  pnm2_latmaxi, pnm2_latmaxi = ADDP(0.0));


                /* Im+2,m, Im+3,m, ..., Inmax,m */
                /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                /* This is not a typo, "pnm1" are indeed required here */
                inm2 = (hm[n] * bnm[n]) * inm0 - (anm[n] / (REAL)(n + 1)) *
                       ((u2 * u2) * pnm1_latmaxi - (u1 * u1) * pnm1_latmini);

                inm0 = inm1;
                inm1 = inm2;
                /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


                /* Lumped coefficients */
                /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
                inm_cnm = inm2 * shcs->c[m][n - m];
                inm_snm = inm2 * shcs->s[m][n - m];


                *a += rpows[n] * inm_cnm;
                *b += rpows[n] * inm_snm;


                if (symmi)
                {
                    if (npm_even)
                    {
                        *a2 += rpows2[n] * inm_cnm;
                        *b2 += rpows2[n] * inm_snm;
                    }
                    else
                    {
                        *a2 -= rpows2[n] * inm_cnm;
                        *b2 -= rpows2[n] * inm_snm;
                    }
                }
                /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */


                pnm1_latmini = pnm2_latmini;
                pnm1_latmaxi = pnm2_latmaxi;


            } /* End of the loop over harmonic degrees */
            /* ............................................................. */


        } /* End of computation of tesseral harmonics */
        /* ----------------------------------------------------------------- */


    } /* End of computation of the lumped coefficients */
    /* --------------------------------------------------------------------- */






    return;
}
