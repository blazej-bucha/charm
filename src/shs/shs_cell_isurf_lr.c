/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _MSC_VER
#   define _USE_MATH_DEFINES
#endif
#include <math.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_cell_isurf_lr)(REAL lon0, REAL dlon, size_t nlon,
                              REAL lc00, REAL lc01,
                              REAL lc10, REAL lc11,
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
        REAL tmp = dlon * lc00;
        for (size_t j = 0; j < nlon; j++)
            fi[j] += tmp;
    }
    else if ((m1 == 0) || (m3 == 0))
    {
        REAL lon1 = lon0 + dlon;
        REAL lon2 = lon1 + dlon;
        REAL k    = (m1 == 0) ? (REAL)m3 : (REAL)m1;
        REAL k_2  = k / ADDP(2.0);
        REAL d    = ADDP(2.0) * COS(k * dlon);


        REAL z  = ADDP(2.0) * SIN(k_2 * dlon) / k;
        REAL lon_tmp = k_2 * (lon0 + lon1);
        REAL zs = SIN(lon_tmp);
        REAL zc = COS(lon_tmp);
        REAL tmp0;
        if (m1 == 0)
            tmp0 = z * (lc00 * zc + lc01 * zs);
        else
            tmp0 = z * (lc00 * zc + lc10 * zs);
        fi[0] += tmp0;


        if (nlon == 1)
            return;


        lon_tmp = k_2 * (lon1 + lon2);
        zs = SIN(lon_tmp);
        zc = COS(lon_tmp);
        REAL tmp1;
        if (m1 == 0)
            tmp1 = z * (lc00 * zc + lc01 * zs);
        else
            tmp1 = z * (lc00 * zc + lc10 * zs);
        fi[1] += tmp1;


        if (nlon == 2)
            return;


        REAL tmp2;
        for (size_t j = 2; j < nlon; j++)
        {
            tmp2 = d * tmp1 - tmp0;


            fi[j] += tmp2;


            tmp0 = tmp1;
            tmp1 = tmp2;
        }
    }
    else if (m1 == m3)
    {
        REAL lon1 = lon0 + dlon;
        REAL lon2 = lon1 + dlon;
        /* Note that "m1 == m3" */
        REAL k    = ADDP(2.0) * (REAL)m1;
        REAL k_2  = k / ADDP(2.0);
        REAL d    = ADDP(2.0) * COS(k * dlon);


        REAL z  = SIN(k_2 * dlon) / k;
        REAL lon_tmp = k_2 * (lon0 + lon1);
        REAL zs = SIN(lon_tmp);
        REAL zc = COS(lon_tmp);
        REAL tmp0 = z * (zs * (lc01 + lc10) + zc * (lc00 - lc11));
        REAL c    = dlon / ADDP(2.0) * (lc00 + lc11);
        fi[0] += c + tmp0;


        if (nlon == 1)
            return;


        lon_tmp = k_2 * (lon1 + lon2);
        zs = SIN(lon_tmp);
        zc = COS(lon_tmp);
        REAL tmp1 = z * (zs * (lc01 + lc10) + zc * (lc00 - lc11));
        fi[1] += c + tmp1;


        if (nlon == 2)
            return;


        REAL tmp2;
        for (size_t j = 2; j < nlon; j++)
        {
            tmp2 = d * tmp1 - tmp0;


            fi[j] += c + tmp2;


            tmp0 = tmp1;
            tmp1 = tmp2;
        }
    }
    else
    {
        REAL lon1 = lon0 + dlon;
        REAL lon2 = lon1 + dlon;
        REAL k1   = (REAL)m3 - (REAL)m1;
        REAL k2   = (REAL)m3 + (REAL)m1;
        REAL k1_2 = k1 / ADDP(2.0);
        REAL k2_2 = k2 / ADDP(2.0);
        REAL d1   = ADDP(2.0) * COS(k1 * dlon);
        REAL d2   = ADDP(2.0) * COS(k2 * dlon);


        REAL z1 = SIN(k1_2 * dlon) / k1;
        REAL z2 = SIN(k2_2 * dlon) / k2;

        REAL lon_tmp = lon0 + lon1;
        REAL lon_tmp_k1_2 = k1_2 * lon_tmp;
        REAL lon_tmp_k2_2 = k2_2 * lon_tmp;
        REAL tmp0_1 = z1 * ((lc00 + lc11) * COS(lon_tmp_k1_2) +
                              (lc01 - lc10) * SIN(lon_tmp_k1_2));
        REAL tmp0_2 = z2 * ((lc00 - lc11) * COS(lon_tmp_k2_2) +
                              (lc01 + lc10) * SIN(lon_tmp_k2_2));
        fi[0] += tmp0_1 + tmp0_2;


        if (nlon == 1)
            return;


        lon_tmp = lon1 + lon2;
        lon_tmp_k1_2 = k1_2 * lon_tmp;
        lon_tmp_k2_2 = k2_2 * lon_tmp;
        REAL tmp1_1 = z1 * ((lc00 + lc11) * COS(lon_tmp_k1_2) +
                              (lc01 - lc10) * SIN(lon_tmp_k1_2));
        REAL tmp1_2 = z2 * ((lc00 - lc11) * COS(lon_tmp_k2_2) +
                              (lc01 + lc10) * SIN(lon_tmp_k2_2));
        fi[1] += tmp1_1 + tmp1_2;


        if (nlon == 2)
            return;


        REAL tmp2_1, tmp2_2;
        for (size_t j = 2; j < nlon; j++)
        {
            tmp2_1 = d1 * tmp1_1 - tmp0_1;
            tmp2_2 = d2 * tmp1_2 - tmp0_2;


            fi[j] += tmp2_1 + tmp2_2;


            tmp0_1 = tmp1_1;
            tmp1_1 = tmp2_1;
            tmp0_2 = tmp1_2;
            tmp1_2 = tmp2_2;
        }
    }






    return;
}
