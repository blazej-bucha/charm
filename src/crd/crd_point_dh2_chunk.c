/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_propagate.h"
#include "crd_point_dh_chunk.h"
#include "crd_point_dh2_chunk.h"
/* ------------------------------------------------------------------------- */







CHARM(point) *CHARM(crd_point_dh2_chunk)(unsigned long nmax,
                                         REAL r,
                                         size_t local_nlat,
                                         size_t local_0_start,
                                         CHARM(err) *err)
{
    CHARM(point) *dhg = CHARM(crd_point_dh_chunk)(nmax, r,
                                                  CHARM_CRD_POINT_GRID_DH2,
                                                  local_nlat, local_0_start,
                                                  err);
    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    return dhg;
}
