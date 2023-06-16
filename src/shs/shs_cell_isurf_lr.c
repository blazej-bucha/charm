/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_isurf_lr)(REAL lon0, REAL dlon, size_t nlon,
                              REAL_SIMD lc00, REAL_SIMD lc01,
                              REAL_SIMD lc10, REAL_SIMD lc11,
                              unsigned long m1,
                              unsigned long m3,
                              REAL *fi)
/*
 * ============================================================================
 *
 * DESCRIPTION: Computes partial sums for the "i"th grid latitude inside the
 *              "CHARM(shs_cell_isurf)" function.
 *
 *
 * INPUTS: "lon0", "dlon", ..., "m3" -- All variables have the same meaning as
 *                                      in "CHARM(shs_cell_isurf)",
 *                                      where further details can be found.
 *
 *
 * OUTPUTS: "fi" -- A pointer to an array of dimensions "nlon" to store the
 *                  partial contributions to the the "i"th grid latitude.
 *
 * ============================================================================
 *
 * */
{
    if ((m1 == 0) && (m3 == 0))
    {
        REAL_SIMD tmp = MUL_R(SET1_R(dlon), lc00);
        for (size_t j = 0; j < nlon; j++)
            STORE_R(&fi[j * SIMD_SIZE], ADD_R(LOAD_R(&fi[j * SIMD_SIZE]),
                                              tmp));
    }
    else if ((m1 == 0) || (m3 == 0))
    {
        REAL lon1   = lon0 + dlon;
        REAL lon2   = lon1 + dlon;
        REAL k      = (m1 == 0) ? (REAL)m3 : (REAL)m1;
        REAL k_2    = k / PREC(2.0);
        REAL_SIMD d = SET1_R(PREC(2.0) * COS(k * dlon));


        REAL_SIMD z  = SET1_R(PREC(2.0) * SIN(k_2 * dlon) / k);
        REAL lon_tmp = k_2 * (lon0 + lon1);
        REAL_SIMD zs = SET1_R(SIN(lon_tmp));
        REAL_SIMD zc = SET1_R(COS(lon_tmp));
        REAL_SIMD tmp0;
        if (m1 == 0)
            tmp0 = MUL_R(z, ADD_R(MUL_R(lc00, zc), MUL_R(lc01, zs)));
        else
            tmp0 = MUL_R(z, ADD_R(MUL_R(lc00, zc), MUL_R(lc10, zs)));
        STORE_R(&fi[0], ADD_R(LOAD_R(&fi[0]), tmp0));


        if (nlon == 1)
            return;


        lon_tmp = k_2 * (lon1 + lon2);
        zs = SET1_R(SIN(lon_tmp));
        zc = SET1_R(COS(lon_tmp));
        REAL_SIMD tmp1;
        if (m1 == 0)
            tmp1 = MUL_R(z, ADD_R(MUL_R(lc00, zc), MUL_R(lc01, zs)));
        else
            tmp1 = MUL_R(z, ADD_R(MUL_R(lc00, zc), MUL_R(lc10, zs)));
        STORE_R(&fi[SIMD_SIZE], ADD_R(LOAD_R(&fi[SIMD_SIZE]), tmp1));


        if (nlon == 2)
            return;


        REAL_SIMD tmp2;
        for (size_t j = 2; j < nlon; j++)
        {
            tmp2 = SUB_R(MUL_R(d, tmp1), tmp0);


            STORE_R(&fi[j * SIMD_SIZE], ADD_R(LOAD_R(&fi[j * SIMD_SIZE]),
                                              tmp2));


            tmp0 = tmp1;
            tmp1 = tmp2;
        }
    }
    else if (m1 == m3)
    {
        REAL lon1   = lon0 + dlon;
        REAL lon2   = lon1 + dlon;
        /* Note that "m1 == m3" */
        REAL k      = PREC(2.0) * (REAL)m1;
        REAL k_2    = k / PREC(2.0);
        REAL_SIMD d = SET1_R(PREC(2.0) * COS(k * dlon));


        REAL_SIMD z    = SET1_R(SIN(k_2 * dlon) / k);
        REAL lon_tmp   = k_2 * (lon0 + lon1);
        REAL_SIMD zs   = SET1_R(SIN(lon_tmp));
        REAL_SIMD zc   = SET1_R(COS(lon_tmp));
        REAL_SIMD tmp0 = MUL_R(z, ADD_R(MUL_R(zs, ADD_R(lc01, lc10)),
                                        MUL_R(zc, SUB_R(lc00, lc11))));
        REAL_SIMD c    = MUL_R(SET1_R(dlon / PREC(2.0)), ADD_R(lc00, lc11));
        STORE_R(&fi[0], ADD_R(LOAD_R(&fi[0]), ADD_R(c, tmp0)));


        if (nlon == 1)
            return;


        lon_tmp = k_2 * (lon1 + lon2);
        zs = SET1_R(SIN(lon_tmp));
        zc = SET1_R(COS(lon_tmp));
        REAL_SIMD tmp1 = MUL_R(z, ADD_R(MUL_R(zs, ADD_R(lc01, lc10)),
                                        MUL_R(zc, SUB_R(lc00, lc11))));
        STORE_R(&fi[SIMD_SIZE], ADD_R(LOAD_R(&fi[SIMD_SIZE]), ADD_R(c, tmp1)));


        if (nlon == 2)
            return;


        REAL_SIMD tmp2;
        for (size_t j = 2; j < nlon; j++)
        {
            tmp2 = SUB_R(MUL_R(d, tmp1), tmp0);


            STORE_R(&fi[j * SIMD_SIZE], ADD_R(LOAD_R(&fi[j * SIMD_SIZE]),
                                              ADD_R(c, tmp2)));


            tmp0 = tmp1;
            tmp1 = tmp2;
        }
    }
    else
    {
        REAL lon1    = lon0 + dlon;
        REAL lon2    = lon1 + dlon;
        REAL k1      = (REAL)m3 - (REAL)m1;
        REAL k2      = (REAL)m3 + (REAL)m1;
        REAL k1_2    = k1 / PREC(2.0);
        REAL k2_2    = k2 / PREC(2.0);
        REAL_SIMD d1 = SET1_R(PREC(2.0) * COS(k1 * dlon));
        REAL_SIMD d2 = SET1_R(PREC(2.0) * COS(k2 * dlon));


        REAL_SIMD z1 = SET1_R(SIN(k1_2 * dlon) / k1);
        REAL_SIMD z2 = SET1_R(SIN(k2_2 * dlon) / k2);


        REAL lon_tmp = lon0 + lon1;
        REAL lon_tmp_k1_2 = k1_2 * lon_tmp;
        REAL lon_tmp_k2_2 = k2_2 * lon_tmp;
        REAL_SIMD tmp0_1 = MUL_R(z1, ADD_R(MUL_R(ADD_R(lc00, lc11),
                                                 SET1_R(COS(lon_tmp_k1_2))),
                                           MUL_R(SUB_R(lc01, lc10),
                                                 SET1_R(SIN(lon_tmp_k1_2)))));
        REAL_SIMD tmp0_2 = MUL_R(z2, ADD_R(MUL_R(SUB_R(lc00, lc11),
                                                 SET1_R(COS(lon_tmp_k2_2))),
                                           MUL_R(ADD_R(lc01, lc10),
                                                 SET1_R(SIN(lon_tmp_k2_2)))));
        STORE_R(&fi[0], ADD_R(LOAD_R(&fi[0]), ADD_R(tmp0_1, tmp0_2)));


        if (nlon == 1)
            return;


        lon_tmp = lon1 + lon2;
        lon_tmp_k1_2 = k1_2 * lon_tmp;
        lon_tmp_k2_2 = k2_2 * lon_tmp;
        REAL_SIMD tmp1_1 = MUL_R(z1, ADD_R(MUL_R(ADD_R(lc00, lc11),
                                                 SET1_R(COS(lon_tmp_k1_2))),
                                           MUL_R(SUB_R(lc01, lc10),
                                                 SET1_R(SIN(lon_tmp_k1_2)))));
        REAL_SIMD tmp1_2 = MUL_R(z2, ADD_R(MUL_R(SUB_R(lc00, lc11),
                                                 SET1_R(COS(lon_tmp_k2_2))),
                                           MUL_R(ADD_R(lc01, lc10),
                                                 SET1_R(SIN(lon_tmp_k2_2)))));
        STORE_R(&fi[SIMD_SIZE], ADD_R(LOAD_R(&fi[SIMD_SIZE]),
                                      ADD_R(tmp1_1, tmp1_2)));


        if (nlon == 2)
            return;


        REAL_SIMD tmp2_1, tmp2_2;
        for (size_t j = 2; j < nlon; j++)
        {
            tmp2_1 = SUB_R(MUL_R(d1, tmp1_1), tmp0_1);
            tmp2_2 = SUB_R(MUL_R(d2, tmp1_2), tmp0_2);


            STORE_R(&fi[j * SIMD_SIZE], ADD_R(LOAD_R(&fi[j * SIMD_SIZE]),
                                              ADD_R(tmp2_1, tmp2_2)));


            tmp0_1 = tmp1_1;
            tmp1_1 = tmp2_1;
            tmp0_2 = tmp1_2;
            tmp1_2 = tmp2_2;
        }
    }






    return;
}
