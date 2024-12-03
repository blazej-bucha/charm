/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <math.h>
#include "../prec.h"
#include "../crd/crd_point_isEmpty.h"
#include "../crd/crd_point_isGrid.h"
#include "../misc/misc_arr_chck_symm.h"
#include "crd_point_issymm.h"
/* ------------------------------------------------------------------------- */






_Bool CHARM(crd_point_issymm)(const CHARM(point) *pnt)
{
    if (CHARM(crd_point_isEmpty)(pnt))
        return 0;


    if (!CHARM(crd_point_isGrid)(pnt->type))
        return 0;


    size_t local_nlat;
#if HAVE_MPI
    local_nlat = pnt->local_nlat;
#else
    local_nlat = pnt->nlat;
#endif


    /* Now "pnt" has at least one latitude */


    CHARM(err) *err = CHARM(err_init)();
    if (err == NULL)
        return 0;


    int ret = CHARM(misc_arr_chck_symm)(pnt->lat, local_nlat, PREC(0.0),
                                        CHARM(glob_threshold), err);
    if (!CHARM(err_isempty)(err))
        ret = 9999;


    CHARM(err_free)(err);


    return ret == 0;
}
