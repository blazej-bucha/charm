/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






/* An internal function to check whether or not FFT can be applied for
 * grid-wise synthesis. */
void CHARM(shs_grd_fft_lc)(unsigned long m, REAL dlon,
                           REAL_SIMD a, REAL_SIMD b,
                           REAL_SIMD a2, REAL_SIMD b2,
                           _Bool symm, REAL_SIMD symm_simd,
                           int grd_type, REAL *lc_tmp, REAL *lc2_tmp)
{
    REAL_SIMD c = (m == 0) ? SET1_R(PREC(1.0)) : SET1_R(PREC(0.5));
    size_t idx  = m * 2 * SIMD_SIZE;
    size_t idx2 = idx + SIMD_SIZE;


    if (grd_type == CHARM_CRD_CELL_GRID)
    {
        /* Due to the use of the mean values, some additional terms need to be
         * taken into account when compared with the harmonic analysis based on
         * point values (see, e.g., Colombo, 1981) */
        REAL_SIMD cm, sm;
        if (m == 0)
        {
            sm = SET1_R(dlon);
            cm = SET_ZERO_R;
        }
        else
        {
            REAL mr = (REAL)m;
            cm = SET1_R((COS(mr * dlon) - PREC(1.0)) / mr);
            sm = SET1_R(SIN(mr * dlon) / mr);
        }


        STORE_R(&lc_tmp[idx], MUL_R(SUB_R(MUL_R(a, sm), MUL_R(b, cm)), c));
        STORE_R(&lc_tmp[idx2], -MUL_R(ADD_R(MUL_R(a, cm),
                                                 MUL_R(b, sm)), c));


        if (symm)
        {
            STORE_R(&lc2_tmp[idx], MUL_R(symm_simd, MUL_R(SUB_R(MUL_R(a2, sm),
                                                                MUL_R(b2, cm)),
                                                          c)));
            STORE_R(&lc2_tmp[idx2],
                    -MUL_R(symm_simd, MUL_R(ADD_R(MUL_R(a2, cm),
                                                  MUL_R(b2, sm)), c)));
        }


        return;
    }
    else if ((grd_type == CHARM_CRD_POINT_GRID) ||
             (grd_type == CHARM_CRD_POINT_GRID_GL) ||
             (grd_type == CHARM_CRD_POINT_GRID_DH1) ||
             (grd_type == CHARM_CRD_POINT_GRID_DH2))
    {
        /* Let's prepare the complex Fourier coefficients */
        STORE_R(&lc_tmp[idx], MUL_R(a, c));
        STORE_R(&lc_tmp[idx2], -MUL_R(b, c));


        if (symm)
        {
            STORE_R(&lc2_tmp[idx], MUL_R(a2, c));
            STORE_R(&lc2_tmp[idx2], -MUL_R(b2, c));
        }
    }


    return;
}
