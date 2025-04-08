/* This header file is not a part of API. */


#ifndef __GFM_CAP_Q_NORM_DIST_H__
#define __GFM_CAP_Q_NORM_DIST_H__


#include <config.h>
#include <mpfr.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(gfm_cap_q_norm_dist)(const mpfr_t,
                                       const mpfr_t,
                                       const mpfr_t,
                                       mpfr_t);


#ifdef __cplusplus
}
#endif


#endif
