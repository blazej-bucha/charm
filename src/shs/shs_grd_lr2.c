/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../crd/crd_cell_isGrid.h"
#include "../glob/glob_get_shs_block_lat_multiplier.h"
#include "shs_grd_lr2.h"
/* ------------------------------------------------------------------------- */






/* Macros */
/* ------------------------------------------------------------------------- */
#undef LOOP1
#define LOOP1   for (j = 0; j < nlon; j++)                                    \
                      f[row + j] = ctmp * fi[j * size_blk + lssv];


#undef LOOP2
#define LOOP2   for (j = 0; j < nlon; j++)                                    \
                       f[row + j] = ctmp * fi2[j * size_blk + lssv];
/* ------------------------------------------------------------------------- */






void CHARM(shs_grd_lr2)(size_t i,
                        const REAL *latsinv,
                        int grd_type,
                        size_t nlat,
                        size_t nlon,
                        const REAL *symmv,
                        REAL c,
                        const REAL *latminv,
                        const REAL *latmaxv,
                        REAL deltalon,
                        const REAL *fi,
                        const REAL *fi2,
                        REAL *f)
{
#if HAVE_MPI
    const size_t BLOCK_S = CHARM(glob_get_shs_block_lat_multiplier)();
#else
#   define BLOCK_S SIMD_BLOCK_S
#endif


    /* Look into "shs_grd_fft" for the explanation of these variables */
    size_t ipv, j, row, lss, lssv;
    _Bool is_cell_grd = CHARM(crd_cell_isGrid)(grd_type);
    size_t simd_blk = (is_cell_grd) ? 1 : BLOCK_S;
    size_t size_blk = SIMD_SIZE * simd_blk;
    REAL dsigma;
    REAL ctmp = c;


    for (size_t l = 0; l < simd_blk; l++)
    {
        lss = l * SIMD_SIZE;


        for (size_t v = 0; v < SIMD_SIZE; v++)
        {
            lssv = lss + v;


            if (latsinv[lssv] == 0)
                continue;


            if (is_cell_grd)
            {
                dsigma = (SIN(latmaxv[lssv]) - SIN(latminv[lssv])) * deltalon;
                ctmp   = c / dsigma;
            }


            ipv = i + lssv;
            row = ipv * nlon;


            /* If "is_cell_grd" is true (shs with cells), the function is
             * called inside a parallel region, so parallelization of "LOOP1"
             * and "LOOP2" is not possible.  Otherwise (shs with points), the
             * loops are called outside a parallel region, so we parallelize
             * the loop. */
            if (is_cell_grd)
            {
                LOOP1;
            }
            else
            {
#if HAVE_OPENMP
#pragma omp parallel for default(shared) private(j)
#endif
                LOOP1;
            }


            if (symmv[lssv])
            {
                row = (nlat - ipv - 1) * nlon;
                if (is_cell_grd)
                {
                    LOOP2;
                }
                else
                {
#if HAVE_OPENMP
#pragma omp parallel for default(shared) private(j)
#endif
                    LOOP2;
                }
            }
        }
    }


    return;
}
