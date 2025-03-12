/* Header files */
/* ------------------------------------------------------------------------- */
#include <config.h>
#include <stdio.h>
#include <mpfr.h>
#include "../prec.h"
#include "../err/err_propagate.h"
#include "../mpfr/mpfr_flush_unreleased_memory.h"
#include "gfm_cap_q_check_type.h"
#include "gfm_cap_qu0.h"
#include "gfm_cap_quu.h"
/* ------------------------------------------------------------------------- */






void CHARM(gfm_cap_q)(const mpfr_t rref,
                      const mpfr_t r,
                      const mpfr_t psi,
                      unsigned long nmax,
                      unsigned pmax,
                      unsigned kmin,
                      unsigned kmax,
                      unsigned imax,
                      int zone,
                      unsigned type,
                      mpfr_prec_t NBITS,
                      mpfr_t *qkpin,
                      CHARM(err) *err)
{
    CHARM(gfm_cap_q_check_type)(type, err);
    if (!CHARM(err_isempty)(err))
    {
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);
        return;
    }


    if (type == CHARM_GFM_Q00)
        CHARM(gfm_cap_qu0)(rref, r, psi, nmax, pmax, kmin, kmax, imax, zone,
                           NBITS, qkpin, err);
    else if (type == CHARM_GFM_Q10)
        CHARM(gfm_cap_qu0)(rref, r, psi, nmax, pmax, kmin + 1, kmax + 1, imax,
                           zone, NBITS, qkpin, err);
    else if (type == CHARM_GFM_Q11)
        CHARM(gfm_cap_quu)(rref, r, psi, nmax, pmax, kmin, kmax, imax, 1, zone,
                           NBITS, qkpin, err);
    else if (type == CHARM_GFM_Q20)
        CHARM(gfm_cap_qu0)(rref, r, psi, nmax, pmax, kmin + 2, kmax + 2, imax,
                           zone, NBITS, qkpin, err);
    else if (type == CHARM_GFM_Q21)
        CHARM(gfm_cap_quu)(rref, r, psi, nmax, pmax, kmin + 1, kmax + 1, imax,
                           1, zone, NBITS, qkpin, err);
    else if (type == CHARM_GFM_Q22)
        CHARM(gfm_cap_quu)(rref, r, psi, nmax, pmax, kmin, kmax, imax, 2, zone,
                           NBITS, qkpin, err);


    if (!CHARM(err_isempty)(err))
        CHARM(err_propagate)(err, __FILE__, __LINE__, __func__);


    mpfr_free_cache();
    FLUSH_UNRELEASED_MEMORY;


    return;
}
