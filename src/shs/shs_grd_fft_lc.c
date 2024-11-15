/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../crd/crd_point_isGrid.h"
#include "../crd/crd_cell_isGrid.h"
#include "../glob/glob_get_shs_block_lat_multiplier.h"
#include "shs_lc_struct.h"
#include "shs_max_npar.h"
#include "shs_point_gradn.h"
#include "shs_grd_fft_lc.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
#undef FC
#define FC(pnfc, par)                                                         \
        for (l = 0; l < simd_blk; l++)                                        \
        {                                                                     \
            STORE_R(&fc_tmp[pnfc + idx + l * SIMD_SIZE],                      \
                    MUL_R(CAT2(lc->a, par, )[l], c));                         \
            STORE_R(&fc_tmp[pnfc + idx2 + l * SIMD_SIZE],                     \
                    NEG_R(MUL_R(CAT2(lc->b, par, )[l], c)));                  \
        }                                                                     \
                                                                              \
                                                                              \
        if (symm)                                                             \
        {                                                                     \
            for (l = 0; l < simd_blk; l++)                                    \
            {                                                                 \
                STORE_R(&fc2_tmp[pnfc + idx + l * SIMD_SIZE],                 \
                        MUL_R(CAT2(lc->a, par, 2)[l], c));                    \
                STORE_R(&fc2_tmp[pnfc + idx2 + l * SIMD_SIZE],                \
                        NEG_R(MUL_R(CAT2(lc->b, par, 2)[l], c)));             \
            }                                                                 \
        }
/* ------------------------------------------------------------------------- */






/* An internal function to prepare the Fourier coefficients for synthesis along
 * the latitude parallels. */
void CHARM(shs_grd_fft_lc)(unsigned long m,
                           REAL deltalon,  /* Step in longitudes */
                           int grad,
                           CHARM(lc) *lc,
                           _Bool symm,
                           REAL_SIMD *symm_simd,
                           int grd_type,
                           size_t nfc,
                           REAL *fc_tmp,
                           REAL *fc2_tmp)
{
#if HAVE_MPI
    const size_t BLOCK_S = CHARM(glob_get_shs_block_lat_multiplier)();
#else
#   define BLOCK_S SIMD_BLOCK_S
#endif


    REAL_SIMD c = (m == 0) ? SET1_R(PREC(1.0)) : SET1_R(PREC(0.5));
    _Bool is_cell_grd = CHARM(crd_cell_isGrid)(grd_type);
    size_t simd_blk = (is_cell_grd) ? 1 : BLOCK_S;
    size_t size_blk = SIMD_SIZE * simd_blk;
    size_t l;
    size_t idx  = m * 2 * size_blk;
    size_t idx2 = idx + size_blk;


    if (is_cell_grd)
    {
        /* Due to the use of the mean values, some additional terms need to be
         * taken into account when compared with point values (see, e.g.,
         * Colombo, 1981) */
        REAL_SIMD cm, sm;
        if (m == 0)
        {
            sm = SET1_R(deltalon);
            cm = SET_ZERO_R;
        }
        else
        {
            REAL mr = (REAL)m;
            cm = SET1_R((COS(mr * deltalon) - PREC(1.0)) / mr);
            sm = SET1_R(SIN(mr * deltalon) / mr);
        }


        /* With cells, "lc->a", "lc->b", "lc->a2" and "lc->b2" implicitly
         * assume that "BLOCK_S" is "1". */
        for (l = 0; l < simd_blk; l++)
        {
            STORE_R(&fc_tmp[idx + l * SIMD_SIZE],
                    MUL_R(SUB_R(MUL_R(lc->a[0], sm), MUL_R(lc->b[0], cm)), c));
            STORE_R(&fc_tmp[idx2 + l * SIMD_SIZE],
                    NEG_R(MUL_R(ADD_R(MUL_R(lc->a[0], cm),
                                      MUL_R(lc->b[0], sm)),
                                c)));
        }


        if (symm)
        {
            for (l = 0; l < simd_blk; l++)
            {
                STORE_R(&fc2_tmp[idx], MUL_R(*symm_simd,
                                             MUL_R(SUB_R(MUL_R(lc->a2[0], sm),
                                                         MUL_R(lc->b2[0], cm)),
                                                   c)));
                STORE_R(&fc2_tmp[idx2],
                        NEG_R(MUL_R(*symm_simd,
                                    MUL_R(ADD_R(MUL_R(lc->a2[0], cm),
                                                MUL_R(lc->b2[0], sm)),
                                          c))));
            }
        }
    }
    else if (CHARM(crd_point_isGrid)(grd_type))
    {
        /* The ordering of the coefficients in "fc_tmp" and "fc2_tmp" is as
         * follows:
         *
         * * if computing only one parameter, the coefficients of these
         *   parameter are first,
         *
         * * if computing "GRAD_1", the ordering is "vl", "vr" and "vp",
         *
         * * if computing "GRAD_2", the ordering is "vll", "vlr", "vlp", "vrr",
         *   "vrp" and "vpp".
         *
         * In this order we call below the "FC" macro ("l", "r", "p" for
         * "GRAD_1"; "ll", "lr", "lp", "rr", "rp" and "pp" for "GRAD_2").  Do
         * not change this order, as this affects the of the output quantities.
         * This is also the order of the symbolic constants "GRAD_L", "GRAD_R",
         * ... "GRAD_PP" in "shs_point_gradn.h"  */
        FC(0, );


        if (grad > 0)
        {
            size_t q = nfc * SIMD_SIZE * BLOCK_S * 2;
            size_t pnfc = q;


            FC(pnfc, r);
            pnfc += q;
            FC(pnfc, p);


            if (grad > 1)
            {
                pnfc += q;
                FC(pnfc, rr);
                pnfc += q;
                FC(pnfc, rp);
                pnfc += q;
                FC(pnfc, pp);
            }
        }
    }


    return;
}
