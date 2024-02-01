/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../misc/misc_is_nearly_equal.h"
#include "shs_r_eq_rref.h"
/* ------------------------------------------------------------------------- */






/* Returns "true" if all elements of "pnt->r" are equal to "shcs->r".
 * Otherwise, "false" is returned. */
_Bool CHARM(shs_r_eq_rref)(const CHARM(point) *pnt, const CHARM(shc) *shcs)
{
    for (size_t i = 0; i < pnt->nlat; i++)
    {
        if (!CHARM(misc_is_nearly_equal)(pnt->r[i], shcs->r,
                                         CHARM(glob_threshold)))
        {
            return 0;
        }
    }


    return 1;
}
