/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "../simd/calloc_aligned.h"
#include "../simd/free_aligned.h"
#include "../glob/glob_get_shs_block_lat_multiplier.h"
#include "shs_lc_struct.h"
#include "shs_lc_init.h"
/* ------------------------------------------------------------------------- */






/* Internal function to properly set up a "CHARM(lc)" struct.  It *must* be
 * called after the "CHARM(lc)" structure is declared and before it is used for
 * the first time. */
CHARM(lc) *CHARM(shs_lc_init)(void)
{
    CHARM(lc) *x = (CHARM(lc) *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                      1,
                                                      sizeof(CHARM(lc)));
    if (x == NULL)
        return NULL;


#if HAVE_MPI
    const size_t BLOCK_S = CHARM(glob_get_shs_block_lat_multiplier)();
    x->_all = (REAL_SIMD *)CHARM(calloc_aligned)(SIMD_MEMALIGN,
                                                 LC_BLOCKS * BLOCK_S,
                                                 sizeof(REAL_SIMD));
    if (x->_all == NULL)
        goto FAILURE;
#else
#   define BLOCK_S SIMD_BLOCK_S
    /* In this case, "x->_all" is allocated statically.  See also the note in
     * "shs_lc_struct.c". */
#endif


    x->a    = &(x->_all[0 * BLOCK_S]);
    x->b    = &(x->_all[1 * BLOCK_S]);
    x->a2   = &(x->_all[2 * BLOCK_S]);
    x->b2   = &(x->_all[3 * BLOCK_S]);


    x->ar   = &(x->_all[4 * BLOCK_S]);
    x->ap   = &(x->_all[5 * BLOCK_S]);
    x->arr  = &(x->_all[6 * BLOCK_S]);
    x->arp  = &(x->_all[7 * BLOCK_S]);
    x->app  = &(x->_all[8 * BLOCK_S]);


    x->br   = &(x->_all[9  * BLOCK_S]);
    x->bp   = &(x->_all[10 * BLOCK_S]);
    x->brr  = &(x->_all[11 * BLOCK_S]);
    x->brp  = &(x->_all[12 * BLOCK_S]);
    x->bpp  = &(x->_all[13 * BLOCK_S]);


    x->ar2  = &(x->_all[14 * BLOCK_S]);
    x->ap2  = &(x->_all[15 * BLOCK_S]);
    x->arr2 = &(x->_all[16 * BLOCK_S]);
    x->arp2 = &(x->_all[17 * BLOCK_S]);
    x->app2 = &(x->_all[18 * BLOCK_S]);


    x->br2  = &(x->_all[19 * BLOCK_S]);
    x->bp2  = &(x->_all[20 * BLOCK_S]);
    x->brr2 = &(x->_all[21 * BLOCK_S]);
    x->brp2 = &(x->_all[22 * BLOCK_S]);
    x->bpp2 = &(x->_all[23 * BLOCK_S]);


    for (size_t i = 0; i < LC_BLOCKS * BLOCK_S; i++)
        x->_all[i] = SET_ZERO_R;


    x->error = 0;


EXIT:
    return x;


FAILURE:
    CHARM(free_aligned)(x->_all);
    CHARM(free_aligned)(x);
    goto EXIT;
}
