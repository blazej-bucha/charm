/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include "../prec.h"
#include "../err/err_propagate.h"
#include "gfm_check_p.h"
#include "gfm_check_kminkmax.h"
/* ------------------------------------------------------------------------- */






size_t CHARM(gfm_cap_nq)(unsigned long nmax,
                         unsigned pmax,
                         unsigned kmin,
                         unsigned kmax,
                         unsigned imax,
                         CHARM(err) *err)
{
    CHARM(gfm_check_p)(pmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return 0;
    }


    CHARM(gfm_check_kminkmax)(kmin, kmax, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return 0;
    }


    return (kmax - kmin + 1) * pmax * (imax + 1) * (nmax + 1);
}

