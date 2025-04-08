/* This header file is not a part of API. */


#ifndef __GFM_CAP_Q_REF_H__
#define __GFM_CAP_Q_REF_H__


#include <config.h>
#include <mpfr.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(gfm_cap_q_ref)(const mpfr_t,
                                 const mpfr_t,
                                 unsigned long,
                                 unsigned,
                                 unsigned,
                                 unsigned,
                                 unsigned,
                                 unsigned,
                                 mpfr_prec_t,
                                 mpfr_t *,
                                 CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
