/* This header file is not a part of API. */


#ifndef __PARTITION_INTERVAL_H__
#define __PARTITION_INTERVAL_H__


#include <config.h>


#ifdef __cplusplus
extern "C"
{
#endif


extern void partition_interval_ulong(unsigned long,
                                     unsigned long,
                                     unsigned long,
                                     unsigned long,
                                     unsigned long *,
                                     unsigned long *);


extern void partition_interval_size_t(size_t ,
                                      size_t ,
                                      size_t ,
                                      size_t ,
                                      size_t  *,
                                      size_t  *);


#ifdef __cplusplus
}
#endif


#endif
