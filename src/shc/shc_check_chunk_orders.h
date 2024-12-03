/* This header file is not a part of API. */


#ifndef __SHC_CHECK_CHUNK_ORDERS_H__
#define __SHC_CHECK_CHUNK_ORDERS_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern int CHARM(shc_check_chunk_orders)(unsigned long,
                                         size_t,
                                         const unsigned long *,
                                         _Bool,
                                         CHARM(err) *err);


#ifdef __cplusplus
}
#endif


#endif
