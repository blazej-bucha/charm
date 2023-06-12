/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <fftw3.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






/* An internal function to perform FFT-based synthesis along the latitude grid
 * parallels.
 *
 * This routine could potentially be improved by using strides in FFTW. */
void CHARM(shs_grd_fft)(size_t i, int grd_type, size_t nlat, size_t nlon,
                        REAL *latsinv, REAL *latminv, REAL *latmaxv, REAL dlon,
                        FFTW(complex) *lc, FFTW(complex) *lc2, size_t nlc,
                        const REAL *lc_tmp, const REAL *lc2_tmp,
                        REAL mur, const FFTW(plan) plan, const REAL *symmv,
                        REAL *ftmp, REAL *ftmp2, REAL *f)
{
    size_t ipv, row, idx, lss, lssv, lssidx;
    _Bool is_cell_grd = grd_type == CHARM_CRD_CELL_GRID;
    size_t simd_blk = (is_cell_grd) ? 1 : SIMD_BLOCK;
    size_t size_blk = SIMD_SIZE * simd_blk;
    REAL dsigma;
    REAL c = mur;


    for (size_t l = 0; l < simd_blk; l++)
    {
        lss = l * SIMD_SIZE;


        for (size_t v = 0; v < SIMD_SIZE; v++)
        {
            /* This is the local index relative to "i" of the latitude parellel
             * that we are processing. */
            lssv = lss + v;


            if (latsinv[lssv] == 0)
                continue;


            /* This is the global index of the latitude parellel that we are
             * processing */
            ipv = i + lssv;


            /* Index of the first longitude in "f" for a given latitude
             * parallel */
            row = ipv * nlon;


            /* If doing the cell synthesis, we have to take into acount not the
             * size of the cell on the unit sphere and divide by it the "mu
             * / r" factor. */
            if (is_cell_grd)
            {
                /* Cell area on the unit sphere and some useful constants */
                dsigma = (SIN(latmaxv[lssv]) - SIN(latminv[lssv])) * dlon;
                c = mur / dsigma;
            }


            /* Prepare the array with the lumped coefficients that will enter
             * FFTW.  We need to take proper values from "lc_tmp", because the
             * coefficients are in "lc_tmp" in blocks of "SIMD_SIZE" and
             * "simd_blk".  I think the FFTW guru API could be used to avoid
             * this code, but it should not be performance critical. */
            for (size_t j = 0; j < nlc; j++)
            {
                idx      = j * 2 * size_blk + v;
                lssidx   = lss + idx;
                lc[j][0] = lc_tmp[lssidx];
                lc[j][1] = lc_tmp[lssidx + size_blk];
            }


            /* Now the the FFT and the use the scale factor of "mu / r" or so
             * */
            FFTW(execute_dft_c2r)(plan, lc, ftmp);
            for (size_t j = 0; j < nlon; j++)
                f[row + j] = c * ftmp[j];


            if (symmv[lssv])
            {
                for (size_t j = 0; j < nlc; j++)
                {
                    idx       = j * 2 * size_blk + v;
                    lssidx    = lss + idx;
                    lc2[j][0] = lc2_tmp[lssidx];
                    lc2[j][1] = lc2_tmp[lssidx + size_blk];
                }


                row = (nlat - ipv - 1) * nlon;
                FFTW(execute_dft_c2r)(plan, lc2, ftmp2);
                for (size_t j = 0; j < nlon; j++)
                    f[row + j] = c * ftmp2[j];
            }
        }
    }


    return;
}
