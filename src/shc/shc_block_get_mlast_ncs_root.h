/* This header file is not a part of API. */


#ifndef __SHC_BLOCK_GET_MLAST_NCS_ROOT_H__
#define __SHC_BLOCK_GET_MLAST_NCS_ROOT_H__


#include <config.h>
#include "../prec.h"
#include "shc_block_struct.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shc_block_get_mlast_ncs_root)(const CHARM(shc) *,
                                                CHARM(shc_block) *,
                                                unsigned long,
                                                unsigned long *,
                                                size_t *,
                                                int *,
                                                CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
