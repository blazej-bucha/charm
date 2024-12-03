/* This header file is not a part of API. */


#ifndef __SHC_INIT_CHUNK_H__
#define __SHC_INIT_CHUNK_H__


#include <config.h>
#include "../prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern CHARM(shc) *CHARM(shc_init_chunk)(unsigned long,
                                         REAL,
                                         REAL,
                                         REAL *,
                                         REAL *,
                                         size_t,
                                         const unsigned long *,
                                         CHARM(err) *);


#ifdef __cplusplus
}
#endif


#endif
