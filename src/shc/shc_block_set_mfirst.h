/* This header file is not a part of API. */


#ifndef __SHC_BLOCK_SET_MFIRST_H__
#define __SHC_BLOCK_SET_MFIRST_H__


#include <config.h>
#include "../prec.h"
#include "shc_block_struct.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern void CHARM(shc_block_set_mfirst)(CHARM(shc_block) *,
                                        const CHARM(shc) *,
                                        unsigned long,
                                        CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
