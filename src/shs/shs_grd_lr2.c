/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_grd_lr2)(size_t i, const REAL *latsinv,
                        int grd_type, size_t nlat, size_t nlon,
                        const REAL *symmv, REAL c,
                        const REAL *latminv, const REAL *latmaxv,
                        REAL dlon, const REAL *fi, const REAL *fi2, REAL *f)
{
    /* Look into "shs_grd_fft" for the explanation of these variables */
    size_t ipv, row, lss, lssv;
    _Bool is_cell_grd = grd_type == CHARM_CRD_CELL_GRID;
    size_t simd_blk = (is_cell_grd) ? 1 : SIMD_BLOCK;
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
                dsigma = (SIN(latmaxv[lssv]) - SIN(latminv[lssv])) * dlon;
                ctmp   = c / dsigma;
            }


            ipv = i + lssv;
            row = ipv * nlon;


            for (size_t j = 0; j < nlon; j++)
                f[row + j] = ctmp * fi[j * size_blk + lssv];


            if (symmv[lssv])
            {
                row = (nlat - ipv - 1) * nlon;
                for (size_t j = 0; j < nlon; j++)
                    f[row + j] = ctmp * fi2[j * size_blk + lssv];
            }
        }
    }


    return;
}
