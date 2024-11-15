/* This header file is not a part of API. */


#ifndef __SHC_BLOCK_SET_COEFFS_H__
#define __SHC_BLOCK_SET_COEFFS_H__


#include <config.h>
#include "../prec.h"
#include "shc_block_struct.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shc_block_set_coeffs)(CHARM(shc) *
#if HAVE_MPI
                                        , CHARM(shc_block) *,
                                        unsigned long,
                                        unsigned long,
                                        CHARM(err) *
#endif
                                       );


#ifdef __cplusplus
}
#endif


#endif
