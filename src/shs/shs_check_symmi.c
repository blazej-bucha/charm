/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* An internal function to determine whether or not the symmetry property of
 * Legendre functions has to be applied for the "i"th latitude. */
_Bool CHARM(shs_check_symmi)(int grd_type, size_t nlatdo, _Bool symm,
                             _Bool even, size_t i)
{
    if (symm == 0)

        /* The grid was already classified as non-symmetric */
        return 0;

    else if (symm == 1 && even == 0 && i == (nlatdo - 1))

        /* For a symmetric grid containing the equator, this ensures that the
         * points on the equator will not be processed twice */
        return 0;

    else if ((grd_type == CHARM_CRD_POINTS_GRID_DH1 ||
              grd_type == CHARM_CRD_POINTS_GRID_DH2) && (i == 0))

        /* For the Driscoll--Healy grids, do not apply the symmetry property at
         * the north pole */
        return 0;

    else

        /* Exploit the symmetry property */
        return 1;

}
