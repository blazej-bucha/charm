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
 * synthesis with point values */
void CHARM(shs_point_kernel)(unsigned long nmax, unsigned long m,
                             const CHARM(shc) *shcs,
                             const REAL *anm, const REAL *bnm,
                             REAL t, const REAL *ps, const int *ips,
                             const REAL *rpows, const REAL *rpows2,
                             _Bool symmi,
                             REAL *a, REAL *b,
                             REAL *a2, REAL *b2)
{
    REAL w, x,y, z;
    int ixy, ix, iy, iz;


    REAL pnm0, pnm1, pnm2;
    REAL pnm_cnm, pnm_snm;


    _Bool npm_even; /* True if "n + m" is even */


    /* Summations over harmonic degree "n" */
    if (m == 0)
    {
        /* Zonal harmonics */
        /* ----------------------------------------------------- */
        /* P00 */
        pnm0    = ADDP(1.0);
        pnm_cnm = pnm0 * shcs->c[0][0];
        *a      = rpows[0] * pnm_cnm;
        *b      = ADDP(0.0);


        if (symmi)
        {
            *a2 = rpows2[0] * pnm_cnm;
            *b2 = ADDP(0.0);
        }


        /* P10 */
        if (nmax >= 1)
        {
            pnm1    = ROOT3 * t;
            pnm_cnm = pnm1 * shcs->c[0][1];
            *a     += rpows[1] * pnm_cnm;


            if (symmi)
                *a2 -= rpows2[1] * pnm_cnm;
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
                pnm2    = anm[n] * t * pnm1 - bnm[n] * pnm0;
                pnm_cnm = pnm2 * shcs->c[0][n];
                *a     += rpows[n] * pnm_cnm;


                if (symmi)
                {
                    if (npm_even)
                        *a2 += rpows2[n] * pnm_cnm;
                    else
                        *a2 -= rpows2[n] * pnm_cnm;
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
        PNM_SECTORIAL_XNUM(x, ix, ps[m - 1], ips[m - 1], pnm0);


        pnm_cnm = pnm0 * shcs->c[m][0];
        pnm_snm = pnm0 * shcs->s[m][0];


        *a = rpows[m] * pnm_cnm;
        *b = rpows[m] * pnm_snm;


        if (symmi)
        {
            *a2 = rpows2[m] * pnm_cnm;
            *b2 = rpows2[m] * pnm_snm;
        }
        /* ----------------------------------------------------- */


        /* Tesseral harmonics */
        /* ----------------------------------------------------- */
        if (m < nmax)
        {
            PNM_SEMISECTORIAL_XNUM(x, y, ix, iy, w, t, anm[m + 1], pnm1);


            pnm_cnm = pnm1 * shcs->c[m][1];
            pnm_snm = pnm1 * shcs->s[m][1];


            *a += rpows[m + 1] * pnm_cnm;
            *b += rpows[m + 1] * pnm_snm;


            if (symmi)
            {
                *a2 -= rpows2[m + 1] * pnm_cnm;
                *b2 -= rpows2[m + 1] * pnm_snm;
            }


            /* Loop over degrees */
            /* ------------------------------------------------- */
            /* Is "n + m" even?  Since we start the loop with "n = m + 2", then
             * the parity of the first "m + 2 + m" is always even.  Then, it
             * changes with every loop iteration. */
            npm_even = 1;
            for (unsigned long n = (m + 2); n <= nmax;
                 n++, npm_even = !npm_even)
            {
                /* Compute tesseral Legendre function */
                PNM_TESSERAL_XNUM(x, y, z, ix, iy, iz,
                                  ixy, w, t, anm[n], bnm[n],
                                  pnm2, continue);


                pnm_cnm = pnm2 * shcs->c[m][n - m];
                pnm_snm = pnm2 * shcs->s[m][n - m];


                *a += rpows[n] * pnm_cnm;
                *b += rpows[n] * pnm_snm;


                if (symmi)
                {
                    if (npm_even)
                    {
                        *a2 += rpows2[n] * pnm_cnm;
                        *b2 += rpows2[n] * pnm_snm;
                    }
                    else
                    {
                        *a2 -= rpows2[n] * pnm_cnm;
                        *b2 -= rpows2[n] * pnm_snm;
                    }
                }


            } /* End of the loop over harmonic degrees */
            /* ------------------------------------------------- */


        } /* End of computation of tesseral harmonics */
        /* ----------------------------------------------------- */


    } /* End of the summation over harmonic degree "n" */
    /* --------------------------------------------------------- */






    return;
}
