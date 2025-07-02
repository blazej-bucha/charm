/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "shs_rpows.h"
/* ------------------------------------------------------------------------- */






/* An internal function to compute integer powers of "(R / r)". */
void CHARM(shs_rpows)(REAL_SIMD r,
                      REAL_SIMD rref,
                      size_t imax,
                      size_t offset,
                      REAL_SIMD *rpows)
{
    REAL_SIMD ratio = DIV_R(rref, r);


    size_t idx0 = 0;
    size_t idx1 = offset;
    rpows[idx0] = SET1_R(PREC(1.0));
    for (size_t i = 1; i <= imax; i++)
    {
        rpows[idx1] = MUL_R(rpows[idx0], ratio);
        idx0 += offset;
        idx1 += offset;
    }


    return;
}
