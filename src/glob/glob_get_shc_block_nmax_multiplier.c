/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "glob_get_shc_block_nmax_multiplier.h"
/* ------------------------------------------------------------------------- */






unsigned long CHARM(glob_get_shc_block_nmax_multiplier)(void)
{
#if HAVE_MPI
    unsigned long ret = CHARM(glob_shc_block_nmax_multiplier);
#else
    unsigned long ret = 1;
#endif
    return (ret > 0) ? ret : 1;
}
