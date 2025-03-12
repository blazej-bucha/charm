/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_set.h"
#include "gfm_check_kminkmax.h"
/* ------------------------------------------------------------------------- */






/* Checks the value of integer topography power. */
void CHARM(gfm_check_kminkmax)(unsigned kmin, unsigned kmax, CHARM(err) *err)
{
    if (kmin > kmax)
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Order of the minimum radial derivative cannot be "
                       "larger than the order of the maximum radial "
                       "derivative.");


    return;
}
