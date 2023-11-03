/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../crd/crd_point_isGrid.h"
#include "../crd/crd_cell_isGrid.h"
/* ------------------------------------------------------------------------- */






/* An internal function to check whether or not FFT can be applied for
 * grid-wise synthesis. */
void CHARM(shs_grd_fft_lc)(unsigned long m, REAL dlon,
                           REAL_SIMD *a, REAL_SIMD *b,
                           REAL_SIMD *a2, REAL_SIMD *b2,
                           _Bool symm, REAL_SIMD *symm_simd,
                           int grd_type, REAL *lc_tmp, REAL *lc2_tmp)
{
    REAL_SIMD c = (m == 0) ? SET1_R(PREC(1.0)) : SET1_R(PREC(0.5));
    _Bool is_cell_grd = CHARM(crd_cell_isGrid)(grd_type);
    size_t simd_blk = (is_cell_grd) ? 1 : SIMD_BLOCK;
    size_t size_blk = SIMD_SIZE * simd_blk;
#ifdef SIMD
    NEG_R_INIT;
#endif
    size_t l;
    size_t idx  = m * 2 * size_blk;
    size_t idx2 = idx + size_blk;


    if (is_cell_grd)
    {
        /* Due to the use of the mean values, some additional terms need to be
         * taken into account when compared with the harmonic synthesis based
         * on point values (see, e.g., Colombo, 1981) */
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


        for (l = 0; l < simd_blk; l++)
        {
            STORE_R(&lc_tmp[idx + l * SIMD_SIZE], MUL_R(SUB_R(MUL_R(*a, sm),
                                                              MUL_R(*b, cm)),
                                                        c));
            STORE_R(&lc_tmp[idx2 + l * SIMD_SIZE],
                    NEG_R(MUL_R(ADD_R(MUL_R(*a, cm), MUL_R(*b, sm)), c)));
        }


        if (symm)
        {
            for (l = 0; l < simd_blk; l++)
            {
                STORE_R(&lc2_tmp[idx], MUL_R(*symm_simd,
                                             MUL_R(SUB_R(MUL_R(*a2, sm),
                                                         MUL_R(*b2, cm)),
                                                   c)));
                STORE_R(&lc2_tmp[idx2],
                        NEG_R(MUL_R(*symm_simd, MUL_R(ADD_R(MUL_R(*a2, cm),
                                                            MUL_R(*b2, sm)),
                                                      c))));
            }
        }
    }
    else if (CHARM(crd_point_isGrid)(grd_type))
    {
        /* Let's prepare the complex Fourier coefficients */
        for (l = 0; l < simd_blk; l++)
        {
            STORE_R(&lc_tmp[idx + l * SIMD_SIZE], MUL_R(a[l], c));
            STORE_R(&lc_tmp[idx2 + l * SIMD_SIZE], NEG_R(MUL_R(b[l], c)));
        }


        if (symm)
        {
            for (l = 0; l < simd_blk; l++)
            {
                STORE_R(&lc2_tmp[idx + l * SIMD_SIZE], MUL_R(a2[l], c));
                STORE_R(&lc2_tmp[idx2 + l * SIMD_SIZE],
                        NEG_R(MUL_R(b2[l], c)));
            }
        }
    }


    return;
}
