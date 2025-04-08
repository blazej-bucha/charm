/* This header file is not a part of API. */


#ifndef __GFM_CAP_Q_DNORM_DIST_H__
#define __GFM_CAP_Q_DNORM_DIST_H__


#include <config.h>
#include <mpfr.h>
#include "../prec.h"
#include "../mpfr/mpfr_ndarray.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(gfm_cap_q_dnorm_dist)(mpfr_ndarray *,
                                        const mpfr_t,
                                        const mpfr_t,
                                        const mpfr_t,
                                        const mpfr_ndarray *,
                                        const mpfr_ndarray *,
                                        unsigned,
                                        mpfr_prec_t,
                                        CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
