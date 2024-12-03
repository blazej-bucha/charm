/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include "../prec.h"
#include "../simd/simd.h"
#include "glob_get_shs_block_lat_multiplier.h"
/* ------------------------------------------------------------------------- */






size_t CHARM(glob_get_shs_block_lat_multiplier)(void)
{
#if HAVE_MPI
    size_t ret = CHARM(glob_shs_block_lat_multiplier);
#else
    /* This case is for users only.  In the CHarm's source files, we always use
     * a macro to get this value in order to have fixed-length-size arrays,
     * fixed number of loop runs, etc. */
    size_t ret = SIMD_BLOCK_S;
#endif
    return (ret > 0) ? ret : 1;
}
