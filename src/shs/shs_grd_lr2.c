/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../simd/simd.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_grd_lr2)(size_t i, size_t v, size_t nlat, size_t nlon,
                        const REAL *symmv, REAL c, const REAL *fi,
                        const REAL *fi2, REAL *f)
{
    size_t ipv     = i + v;
    size_t row = ipv * nlon;


    for (size_t j = 0; j < nlon; j++)
        f[row + j] = c * fi[j * SIMD_SIZE + v];


    if (symmv[v])
    {
        row = (nlat - ipv - 1) * nlon;
        for (size_t j = 0; j < nlon; j++)
            f[row + j] = c * fi2[j * SIMD_SIZE + v];
    }


    return;
}
