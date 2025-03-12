/* This header file is not a part of API. */


#ifndef __GFM_CAP_Q_CHECK_RADIUS_H__
#define __GFM_CAP_Q_CHECK_RADIUS_H__


#include <config.h>
#include <mpfr.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(gfm_cap_q_check_radius)(const mpfr_t,
                                          mpfr_prec_t,
                                          CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
