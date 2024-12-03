/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../simd/free_aligned.h"
#include "shs_lc_struct.h"
#include "shs_lc_free.h"
/* ------------------------------------------------------------------------- */






void CHARM(shs_lc_free)(CHARM(lc) *x)
{
    if (x == NULL)
        return;


#if HAVE_MPI
    CHARM(free_aligned)(x->_all);
#endif
    CHARM(free_aligned)(x);


    return;
}
