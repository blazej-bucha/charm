/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include "../prec.h"
#include "integ_cc.h"
#include "integ_cs.h"
#include "integ_ss.h"
#include "integ_sc.h"
#include "../misc/misc_nan.h"
#include "../leg/leg_pnmj_check_ordering.h"
#include "../err/err_set.h"
#include "../err/err_propagate.h"
#include "../err/err_check_distribution.h"
/* ------------------------------------------------------------------------- */






REAL CHARM(integ_yi1n1m1yi2n2m2)(REAL cltmin,
                                 REAL cltmax,
                                 REAL lonmin,
                                 REAL lonmax,
                                 _Bool i1,
                                 unsigned long n1,
                                 unsigned long m1,
                                 _Bool i2,
                                 unsigned long n2,
                                 unsigned long m2,
                                 const CHARM(pnmj) *pnmj,
                                 CHARM(err) *err)
{
    /* Some simple error checks */
    /* --------------------------------------------------------------------- */
    CHARM(err_check_distribution)(err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return NAN;
    }


    if (cltmin > cltmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"cltmin\" cannot be larger than \"cltmax\".");
        return NAN;
    }


    if (lonmin > lonmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"lonmin\" cannot be larger than \"lonmax\".");
        return NAN;
    }


    if (n1 > pnmj->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n1\" cannot be larger than \"pnmj->nmax\".");
        return NAN;
    }


    if (n2 > pnmj->nmax)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"n2\" cannot be larger than \"pnmj->nmax\".");
        return NAN;
    }


    if (m1 > n1)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"m1\" cannot be larger than \"n1\".");
        return NAN;
    }


    if (m2 > n2)
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "\"m2\" cannot be larger than \"n2\".");
        return NAN;
    }


    if (CHARM(leg_pnmj_check_ordering)(pnmj->ordering))
    {
        CHARM(err_set)(err, __FILE__, __LINE__, __func__, CHARM_EFUNCARG,
                       "Unsupported value of \"pnmj->ordering\".");
        return NAN;
    }

    /* --------------------------------------------------------------------- */


    /* Integrals of Legendre functions */
    /* --------------------------------------------------------------------- */
    REAL ip = CHARM(integ_pn1m1pn2m2)(cltmin, cltmax, n1, m1, n2, m2, pnmj,
                                      err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return NAN;
    }
    /* --------------------------------------------------------------------- */


    /* Integrals of the trigonometric functions containing the longitude */
    /* --------------------------------------------------------------------- */
    REAL il;
    if (i1 == 0 && i2 == 0)
        CHARM(integ_cc)(lonmin, lonmax - lonmin, 1, (REAL)m1, (REAL)m2, &il);
    else if (i1 == 0 && i2 == 1)
        CHARM(integ_cs)(lonmin, lonmax - lonmin, 1, (REAL)m1, (REAL)m2, &il);
    else if (i1 == 1 && i2 == 0)
        CHARM(integ_sc)(lonmin, lonmax - lonmin, 1, (REAL)m1, (REAL)m2, &il);
    else /* if (i1 == 1 && i2 == 1) */
        CHARM(integ_ss)(lonmin, lonmax - lonmin, 1, (REAL)m1, (REAL)m2, &il);
    /* --------------------------------------------------------------------- */


    /* The integral of spherical harmonic functions */
    /* --------------------------------------------------------------------- */
    return (ip * il);
    /* --------------------------------------------------------------------- */
}
