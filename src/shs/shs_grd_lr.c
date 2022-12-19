/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






/* An internal function to perform the recursion over grid longitudes of the
 * "i"th latitude parallel (the "LR" part of the "PSLR" algorithm). */
void CHARM(shs_grd_lr)(unsigned long m, REAL lon0, REAL dlon, size_t nlon,
                       int grd_type, REAL_SIMD a, REAL_SIMD b, REAL_SIMD a2,
                       REAL_SIMD b2, _Bool symm, REAL *fi, REAL *fi2)
{
    if (grd_type == CHARM_CRD_CELL_GRID)
    {
        /* Before applying the PSLR algorithm from Balmino et al. (2012), the
         * lumped coefficients need to be multiplied by additional terms, which
         * stem from the integration in the longitudinal direction */


        REAL m2sm2dl;
        if (m == 0)
            m2sm2dl = dlon;
        else
        {
            REAL m2 = PREC(2.0) / (REAL)m;
            m2sm2dl = m2 * SIN(dlon / m2);
        }


        a = MUL_R(a, SET1_R(m2sm2dl));
        b = MUL_R(b, SET1_R(m2sm2dl)); /* Remember that "b = 0.0"
                                        * for "m == 0.0", so it can safely be
                                        * multiplied by "m2sm2dl * = dlon" */

        if (symm)
        {
            a2 = MUL_R(a2, SET1_R(m2sm2dl));
            b2 = MUL_R(b2, SET1_R(m2sm2dl));  /* Remember that "b2
                                               * = 0.0" for "m == 0.0", so it
                                               * can safely be multiplied * by
                                               * "m2sm2dl = dlon" */
        }
    }


    /* The PSLR algorithm from Balmino et al. (2012) */
    /* --------------------------------------------------------------------- */
    /* The first longitude point/cell */
    /* ..................................................................... */
    REAL       lontmp = (REAL)m * lon0;
    REAL_SIMD clontmp = SET1_R(COS(lontmp));
    REAL_SIMD slontmp = SET1_R(SIN(lontmp));
    REAL_SIMD     dm0 = ADD_R(MUL_R(a, clontmp), MUL_R(b, slontmp));
    STORE_R(&fi[0], ADD_R(LOAD_R(&fi[0]), dm0));


    REAL_SIMD dm02;
    if (symm)
    {
        dm02 = ADD_R(MUL_R(a2, clontmp), MUL_R(b2, slontmp));
        STORE_R(&fi2[0], ADD_R(LOAD_R(&fi2[0]), dm02));
    }


    if (nlon == 1)
        return;
    /* ..................................................................... */


    /* The second longitude point/cell */
    /* ..................................................................... */
     lontmp       = (REAL)m * (lon0 + dlon);
    clontmp       = SET1_R(COS(lontmp));
    slontmp       = SET1_R(SIN(lontmp));
    REAL_SIMD dm1 = ADD_R(MUL_R(a, clontmp), MUL_R(b, slontmp));
    STORE_R(&fi[SIMD_SIZE], ADD_R(LOAD_R(&fi[SIMD_SIZE]), dm1));


    REAL_SIMD dm12;
    if (symm)
    {
        dm12 = ADD_R(MUL_R(a2, clontmp), MUL_R(b2, slontmp));
        STORE_R(&fi2[SIMD_SIZE], ADD_R(LOAD_R(&fi2[SIMD_SIZE]), dm12));
    }


    if (nlon == 2)
        return;
    /* ..................................................................... */


    /* The third and all the remaining longitude point/cells */
    /* ..................................................................... */
    REAL_SIMD cmdlon2 = SET1_R(PREC(2.0) * COS((REAL)m * dlon));
    REAL_SIMD dm2;
    for (size_t j = 2; j < nlon; j++)
    {
        dm2 = SUB_R(MUL_R(cmdlon2, dm1), dm0);
        STORE_R(&fi[j * SIMD_SIZE], ADD_R(LOAD_R(&fi[j * SIMD_SIZE]), dm2));
        dm0 = dm1;
        dm1 = dm2;
    }


    REAL_SIMD dm22;
    if (symm)
    {
        for (size_t j = 2; j < nlon; j++)
        {
            dm22 = SUB_R(MUL_R(cmdlon2, dm12), dm02);
            STORE_R(&fi2[j * SIMD_SIZE],
                    ADD_R(LOAD_R(&fi2[j * SIMD_SIZE]), dm22));
            dm02 = dm12;
            dm12 = dm22;
        }
    }


    return;
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */


}
