/* This header file is not a part of API. */


#ifndef __SHC_BLOCK_GET_COEFFS_H__
#define __SHC_BLOCK_GET_COEFFS_H__


#include <config.h>
#include "../prec.h"
#include "shc_block_struct.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shc_block_get_coeffs)(const CHARM(shc) *
#if HAVE_MPI
                                        , CHARM(shc_block) *,
                                        unsigned long,
                                        CHARM(err) *
#endif
                                       );


#ifdef __cplusplus
}
#endif


#endif
