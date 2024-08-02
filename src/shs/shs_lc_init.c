/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "shs_lc_struct.h"
#include "shs_lc_init.h"
/* ------------------------------------------------------------------------- */






/* Internal function to properly set up a "CHARM(lc)" struct.  It *must* be
 * called after the "CHARM(lc)" structure is declared and before it is used for
 * the first time. */
void CHARM(shs_lc_init)(CHARM(lc) *x)
{
    x->a    = &(x->_all[0 * SIMD_BLOCK_S]);
    x->b    = &(x->_all[1 * SIMD_BLOCK_S]);
    x->a2   = &(x->_all[2 * SIMD_BLOCK_S]);
    x->b2   = &(x->_all[3 * SIMD_BLOCK_S]);


    x->ar   = &(x->_all[4 * SIMD_BLOCK_S]);
    x->ap   = &(x->_all[5 * SIMD_BLOCK_S]);
    x->arr  = &(x->_all[6 * SIMD_BLOCK_S]);
    x->arp  = &(x->_all[7 * SIMD_BLOCK_S]);
    x->app  = &(x->_all[8 * SIMD_BLOCK_S]);


    x->br   = &(x->_all[9  * SIMD_BLOCK_S]);
    x->bp   = &(x->_all[10 * SIMD_BLOCK_S]);
    x->brr  = &(x->_all[11 * SIMD_BLOCK_S]);
    x->brp  = &(x->_all[12 * SIMD_BLOCK_S]);
    x->bpp  = &(x->_all[13 * SIMD_BLOCK_S]);


    x->ar2  = &(x->_all[14 * SIMD_BLOCK_S]);
    x->ap2  = &(x->_all[15 * SIMD_BLOCK_S]);
    x->arr2 = &(x->_all[16 * SIMD_BLOCK_S]);
    x->arp2 = &(x->_all[17 * SIMD_BLOCK_S]);
    x->app2 = &(x->_all[18 * SIMD_BLOCK_S]);


    x->br2  = &(x->_all[19 * SIMD_BLOCK_S]);
    x->bp2  = &(x->_all[20 * SIMD_BLOCK_S]);
    x->brr2 = &(x->_all[21 * SIMD_BLOCK_S]);
    x->brp2 = &(x->_all[22 * SIMD_BLOCK_S]);
    x->bpp2 = &(x->_all[23 * SIMD_BLOCK_S]);


    for (size_t i = 0; i < LC_BLOCKS * SIMD_BLOCK_S; i++)
        x->_all[i] = SET_ZERO_R;


    return;
}
