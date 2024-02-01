/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../prec.h"
#include "shs_point_grd.h"
#include "shs_point_sctr.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../crd/crd_point_isSctr.h"
#include "../crd/crd_point_isGrid.h"
#include "shs_point_gradn.h"
/* ------------------------------------------------------------------------- */






/* This source file is included in "shs_point_grads.c" and should not be
 * compiled directly.  To avoid this, we compile it only "#if COMPILE_GRADS" is
 * true. */
#if COMPILE_GRADS


#if GRADN == GRAD_0
void CHARM(shs_point_grad0)
#elif GRADN == GRAD_1
void CHARM(shs_point_grad1)
#elif GRADN == GRAD_2
void CHARM(shs_point_grad2)
#else
#   error "GRADN can only be 0, -1 or -2."
#endif
                           (const CHARM(point) *pnt,
                            const CHARM(shc) *shcs,
                            unsigned long nmax,
                            REAL **f,
                            CHARM(err) *err)
{
    /* Some trivial initial error checks */
    /* --------------------------------------------------------------------- */
    if (nmax > shcs->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Maximum harmonic degree of the synthesis (\"nmax\") "
                       "cannot be larger than maximum harmonic degree of "
                       "spherical harmonic coefficients (\"shcs->nmax\").");
        return;
    }
    /* --------------------------------------------------------------------- */






    /* Now do the synthesis */
    /* --------------------------------------------------------------------- */
    if (CHARM(crd_point_isSctr)(pnt->type))
    {
        if (pnt->nlat != pnt->nlon)
        {
            CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                           "The number of latitudes and longitudes in the "
                           "\"pnt\" structure must be the same to "
                           "perform point-wise spherical harmonic synthesis.");
            return;
        }


        /* Point-wise synthesis */
        CHARM(shs_point_sctr)(pnt, shcs, nmax, GRADN, GRADN, GRADN, f, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
    else if (CHARM(crd_point_isGrid)(pnt->type))
    {
        /* Grid-wise synthesis */
        CHARM(shs_point_grd)(pnt, shcs, nmax, GRADN, GRADN, GRADN, f, err);
        if (!CHARM(err_isempty)(err))
        {
            CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
            return;
        }
    }
    else
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported \"pnt->type\" for spherical harmonic "
                       "synthesis of point values.");
        return;
    }
    /* --------------------------------------------------------------------- */






#if GRADN == GRAD_1
    /* For numerical reasons, the gradient elements in "f" are not ordered
     * nicely.  Therefore, swap now the pointers in "f", so that "fx", "fy" and
     * "fz" elements are stored nicely as "f[0]", "f[1]" and "f[2]". */
    /* --------------------------------------------------------------------- */
    REAL *fx = f[GRAD_P];
    REAL *fy = f[GRAD_L];
    REAL *fz = f[GRAD_R];


    f[0] = fx;
    f[1] = fy;
    f[2] = fz;
    /* --------------------------------------------------------------------- */
#elif GRADN == GRAD_2
    /* Do the same but this time with the tensor elements */
    /* --------------------------------------------------------------------- */
    REAL *fxx = f[GRAD_PP];
    REAL *fxy = f[GRAD_LP];
    REAL *fxz = f[GRAD_RP];
    REAL *fyy = f[GRAD_LL];
    REAL *fyz = f[GRAD_LR];
    REAL *fzz = f[GRAD_RR];


    f[0] = fxx;
    f[1] = fxy;
    f[2] = fxz;
    f[3] = fyy;
    f[4] = fyz;
    f[5] = fzz;
    /* --------------------------------------------------------------------- */
#endif






    return;
}


#endif
