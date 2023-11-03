/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../crd/crd_cell_isGrid.h"
/* ------------------------------------------------------------------------- */






/* An internal function to perform the recursion over grid longitudes of the
 * "i"th latitude parallel (the "LR" part of the "PSLR" algorithm). */
void CHARM(shs_grd_lr)(unsigned long m, REAL lon0, REAL dlon, size_t nlon,
                       int grd_type, REAL_SIMD *a, REAL_SIMD *b, REAL_SIMD *a2,
                       REAL_SIMD *b2, _Bool symm, REAL *fi, REAL *fi2)
{
    _Bool is_cell_grd = CHARM(crd_cell_isGrid)(grd_type);
    size_t simd_blk = is_cell_grd ? 1 : SIMD_BLOCK;
    size_t size_blk = SIMD_SIZE * simd_blk;
    size_t lss, jsize_blk;


    if (is_cell_grd)
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


        *a = MUL_R(*a, SET1_R(m2sm2dl));
        *b = MUL_R(*b, SET1_R(m2sm2dl)); /* Remember that "b = 0.0"
                                          * for "m == 0.0", so it can safely be
                                          * multiplied by "m2sm2dl * = dlon" */

        if (symm)
        {
            *a2 = MUL_R(*a2, SET1_R(m2sm2dl));
            *b2 = MUL_R(*b2, SET1_R(m2sm2dl));  /* Remember that "b2
                                                 * = 0.0" for "m == 0.0", so it
                                                 * can safely be multiplied
                                                 * * by "m2sm2dl = dlon" */
        }
    }


    /* The PSLR algorithm from Balmino et al. (2012) */
    /* --------------------------------------------------------------------- */
    /* The first longitude point/cell */
    /* ..................................................................... */
    REAL       lontmp = (REAL)m * lon0;
    REAL_SIMD clontmp = SET1_R(COS(lontmp));
    REAL_SIMD slontmp = SET1_R(SIN(lontmp));
    REAL_SIMD dm0[simd_blk];
    for (size_t l = 0; l < simd_blk; l++)
    {
        lss = l * SIMD_SIZE;
        dm0[l] = ADD_R(MUL_R(a[l], clontmp), MUL_R(b[l], slontmp));
        STORE_R(&fi[lss], ADD_R(LOAD_R(&fi[lss]), dm0[l]));
    }


    REAL_SIMD dm02[simd_blk];
    if (symm)
    {
        for (size_t l = 0; l < simd_blk; l++)
        {
            lss = l * SIMD_SIZE;
            dm02[l] = ADD_R(MUL_R(a2[l], clontmp), MUL_R(b2[l], slontmp));
            STORE_R(&fi2[lss], ADD_R(LOAD_R(&fi2[lss]), dm02[l]));
        }
    }


    if (nlon == 1)
        return;
    /* ..................................................................... */


    /* The second longitude point/cell */
    /* ..................................................................... */
     lontmp       = (REAL)m * (lon0 + dlon);
    clontmp       = SET1_R(COS(lontmp));
    slontmp       = SET1_R(SIN(lontmp));
    REAL_SIMD dm1[simd_blk];
    for (size_t l = 0; l < simd_blk; l++)
    {
        lss = l * SIMD_SIZE;
        dm1[l] = ADD_R(MUL_R(a[l], clontmp), MUL_R(b[l], slontmp));
        STORE_R(&fi[size_blk + lss],
                ADD_R(LOAD_R(&fi[size_blk + lss]),
                      dm1[l]));
    }


    REAL_SIMD dm12[simd_blk];
    if (symm)
    {
        for (size_t l = 0; l < simd_blk; l++)
        {
            lss = l * SIMD_SIZE;
            dm12[l] = ADD_R(MUL_R(a2[l], clontmp), MUL_R(b2[l], slontmp));
            STORE_R(&fi2[size_blk + lss],
                    ADD_R(LOAD_R(&fi2[size_blk + lss]),
                          dm12[l]));
        }
    }


    if (nlon == 2)
        return;
    /* ..................................................................... */


    /* The third and all the remaining longitude points/cells */
    /* ..................................................................... */
    REAL_SIMD cmdlon2 = SET1_R(PREC(2.0) * COS((REAL)m * dlon));
    REAL_SIMD dm2[simd_blk];
    for (size_t j = 2; j < nlon; j++)
    {
        jsize_blk = j * size_blk;


        for (size_t l = 0; l < simd_blk; l++)
        {
            lss = l * SIMD_SIZE;
            dm2[l] = SUB_R(MUL_R(cmdlon2, dm1[l]), dm0[l]);
            STORE_R(&fi[jsize_blk + lss],
                    ADD_R(LOAD_R(&fi[jsize_blk + lss]),
                          dm2[l]));
            dm0[l] = dm1[l];
            dm1[l] = dm2[l];
        }
    }


    REAL_SIMD dm22[simd_blk];
    if (symm)
    {
        for (size_t j = 2; j < nlon; j++)
        {
            jsize_blk = j * size_blk;


            for (size_t l = 0; l < simd_blk; l++)
            {
                lss = l * SIMD_SIZE;
                dm22[l] = SUB_R(MUL_R(cmdlon2, dm12[l]), dm02[l]);
                    STORE_R(&fi2[jsize_blk + lss],
                            ADD_R(LOAD_R(&fi2[jsize_blk + lss]),
                                  dm22[l]));
                dm02[l] = dm12[l];
                dm12[l] = dm22[l];
            }
        }
    }


    return;
    /* ..................................................................... */
    /* --------------------------------------------------------------------- */


}
