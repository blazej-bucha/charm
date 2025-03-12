/* This header file is not a part of API. */


#ifndef __GFM_CAP_Q_DR_H__
#define __GFM_CAP_Q_DR_H__


#include <config.h>
#include <mpfr.h>
#include "../prec.h"
#include "../mpfr/mpfr_ndarray.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(gfm_cap_q_dr)(mpfr_ndarray *,
                                const mpfr_ndarray *,
                                unsigned,
                                unsigned,
                                mpfr_prec_t,
                                CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
