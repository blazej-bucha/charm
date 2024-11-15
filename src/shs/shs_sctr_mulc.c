/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../crd/crd_cell_isSctr.h"
#include "../glob/glob_get_shs_block_lat_multiplier.h"
#include "shs_sctr_mulc.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_sctr_mulc)(size_t i,
                          size_t n,
                          int type,
                          REAL mur,
                          REAL_SIMD tmp,
                          REAL *tmpv,
                          REAL_SIMD *fi,
                          REAL *f)
{
#if HAVE_MPI
    const size_t BLOCK_S = CHARM(glob_get_shs_block_lat_multiplier)();
#else
#   define BLOCK_S SIMD_BLOCK_S
#endif


    size_t ipv;
    size_t il;
    size_t simd_blk = CHARM(crd_cell_isSctr)(type) ? 1 : BLOCK_S;
    REAL_SIMD c = SET1_R(mur);


    for (size_t l = 0; l < simd_blk; l++)
    {
        il = i + l * SIMD_SIZE;
        tmp = MUL_R(c, fi[l]);
        if ((il + SIMD_SIZE) <= n)
            STOREU_R(&f[il], tmp);
        else
        {
            STORE_R(&tmpv[0], tmp);
            for (size_t v = 0; v < SIMD_SIZE; v++)
            {
                ipv = il + v;
                if (ipv < n)
                    f[ipv] = tmpv[v];
                else
                    continue;
            }
        }
    }


    return;
}
