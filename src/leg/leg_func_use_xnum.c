/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdlib.h>
#include "../prec.h"
/* ------------------------------------------------------------------------- */






/* Returns "1" if all elements of "ds" are "1" and "0" otherwise.
 *
 * This function is helpful to check whether the dynamical switching can be 
 * applied during the computation of tesseral Legendre functions.  If *all* 
 * elements of "ds" are "1", this means that all of the one degree higher 
 * Legendre function(s) within the block of latitudes can safely be computed 
 * with the usual "F"-numbers to gain some speed.  If at least one element is 
 * "0", "F"-numbers cannot be applied for all latitudes within the block, so 
 * returned is "0". */
_Bool CHARM(leg_func_use_xnum)(_Bool *ds, size_t nds)
{
    for (size_t i = 0; i < nds; i++)
        /* If "ds[i]" is "0", "X"-numbers have to be used, so no dynamical 
         * switching */
        if (!ds[i])
            return 1;


    return 0;
}
