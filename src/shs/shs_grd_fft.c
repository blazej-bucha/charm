/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../crd/crd_cell_isGrid.h"
#include "shs_grd_fft.h"
/* ------------------------------------------------------------------------- */






/* An internal function to perform FFT-based synthesis along the latitude grid
 * parallels.
 *
 * This routine could potentially be improved by using strides in FFTW. */
void CHARM(shs_grd_fft)(size_t i,
                        int grd_type,
                        size_t nlat,
                        size_t nlon,
                        REAL *latsinv,
                        REAL *latminv,
                        REAL *latmaxv,
                        REAL deltalon,
                        FFTWC(complex) *fc,
                        FFTWC(complex) *fc2,
                        size_t nfc,
                        const REAL *fc_tmp,
                        const REAL *fc2_tmp,
                        REAL mur,
                        const FFTW(plan) plan,
                        const REAL *symmv,
                        REAL *ftmp,
                        REAL *f)
{
    size_t ipv, row, idx, lss, lssv, lssidx;
    _Bool is_cell_grd = CHARM(crd_cell_isGrid)(grd_type);
    size_t simd_blk = (is_cell_grd) ? 1 : SIMD_BLOCK_S;
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


            /* If doing the cell synthesis, we have to take into account the
             * size of the cell on the unit sphere and divide by it the "mu
             * / r" factor. */
            if (is_cell_grd)
            {
                /* Cell area on the unit sphere and some useful constants */
                dsigma = (SIN(latmaxv[lssv]) - SIN(latminv[lssv])) * deltalon;
                c      = mur / dsigma;
            }


            /* Prepare the array with the lumped coefficients that will enter
             * FFTW.  We need to take proper values from "fc_tmp", because the
             * coefficients are in "fc_tmp" in blocks of "SIMD_SIZE" and
             * "simd_blk".  I think the FFTW guru API could be used to avoid
             * this code, but it should not be performance critical. */
            for (size_t j = 0; j < nfc; j++)
            {
                idx      = j * 2 * size_blk + v;
                lssidx   = lss + idx;
                fc[j][0] = fc_tmp[lssidx];
                fc[j][1] = fc_tmp[lssidx + size_blk];
            }


            /* Now execute the FFT and use the scale factor of "mu / r" or so
             * */
            FFTW(execute_dft_c2r)(plan, fc, ftmp);
            for (size_t j = 0; j < nlon; j++)
                f[row + j] = c * ftmp[j];


            if (symmv[lssv])
            {
                for (size_t j = 0; j < nfc; j++)
                {
                    idx       = j * 2 * size_blk + v;
                    lssidx    = lss + idx;
                    fc2[j][0] = fc2_tmp[lssidx];
                    fc2[j][1] = fc2_tmp[lssidx + size_blk];
                }


                row = (nlat - ipv - 1) * nlon;
                FFTW(execute_dft_c2r)(plan, fc2, ftmp);
                for (size_t j = 0; j < nlon; j++)
                    f[row + j] = c * ftmp[j];
            }
        }
    }


    return;
}
