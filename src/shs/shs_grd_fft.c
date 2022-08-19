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
void CHARM(shs_grd_fft)(size_t i, size_t v, size_t nlat, size_t nlon,
                        FFTW(complex) *lc, FFTW(complex) *lc2, size_t nlc,
                        const REAL *lc_tmp, const REAL *lc2_tmp,
                        REAL c, const FFTW(plan) plan, const REAL *symmv,
                        REAL *ftmp, REAL *ftmp2, REAL *f)
{
    size_t ipv = i + v;
    size_t row = ipv * nlon;
    size_t idx;


    for (size_t j = 0; j < nlc; j++)
    {
        idx = j * 2 * SIMD_SIZE + v;
        lc[j][0] = lc_tmp[idx];
        lc[j][1] = lc_tmp[idx + SIMD_SIZE];
    }


    FFTW(execute_dft_c2r)(plan, lc, ftmp);
    for (size_t j = 0; j < nlon; j++)
        f[row + j] = c * ftmp[j];


    if (symmv[v])
    {
        for (size_t j = 0; j < nlc; j++)
        {
            idx = j * 2 * SIMD_SIZE + v;
            lc2[j][0] = lc2_tmp[idx];
            lc2[j][1] = lc2_tmp[idx + SIMD_SIZE];
        }


        row = (nlat - ipv - 1) * nlon;
        FFTW(execute_dft_c2r)(plan, lc2, ftmp2);
        for (size_t j = 0; j < nlon; j++)
            f[row + j] = c * ftmp2[j];
    }


    return;
}
