/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <string.h>
#include "../prec.h"
#include "shc_reset_coeffs.h"
/* ------------------------------------------------------------------------- */






/* An internal function to set all coefficients of a "shc" structure to
 * zero. */
void CHARM(shc_reset_coeffs)(CHARM(shc) *shcs)
{
    size_t idx, nc, ns;
#if HAVE_MPI
    if (shcs->local_nchunk == 0)
        return;


    nc  = shcs->local_nc;
    ns  = shcs->local_ns;
    idx = shcs->local_order[0];
#else
    idx = 0;
    nc  = shcs->nc;
    ns  = shcs->ns;
#endif


    memset(shcs->c[idx], 0, nc * sizeof(REAL));
    memset(shcs->s[idx], 0, ns * sizeof(REAL));


    return;
}
