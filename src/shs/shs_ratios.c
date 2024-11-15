/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../glob/glob_get_shs_block_lat_multiplier.h"
#include "shs_ratios.h"
/* ------------------------------------------------------------------------- */






/* Computes "(R / r)^(m + 1)" ("m" is not typo) */
void CHARM(shs_ratios)(const REAL_SIMD *ratio,
                       const REAL_SIMD *ratio2,
                       _Bool symm,
                       unsigned long m,
                       REAL_SIMD *ratiom,
                       REAL_SIMD *ratio2m)
{
#if HAVE_MPI
    const size_t BLOCK_S = CHARM(glob_get_shs_block_lat_multiplier)();
#else
#   define BLOCK_S SIMD_BLOCK_S
#endif


    size_t l;


    for (l = 0; l < BLOCK_S; l++)
        ratiom[l] = ratio[l];
    if (symm)
        for (l = 0; l < BLOCK_S; l++)
            ratio2m[l] = ratio2[l];


    for (unsigned long mtmp = 1; mtmp <= m; mtmp++)
        for (l = 0; l < BLOCK_S; l++)
            ratiom[l] = MUL_R(ratiom[l], ratio[l]);
    if (symm)
        for (unsigned long mtmp = 1; mtmp <= m; mtmp++)
            for (l = 0; l < BLOCK_S; l++)
                ratio2m[l] = MUL_R(ratio2m[l], ratio2[l]);


    return;
}
