/* This header file is not a part of API. */


#ifndef __WRITE_ARRAY_H__
#define __WRITE_ARRAY_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int write_array(REAL *,
                            size_t,
                            FILE *);


#ifdef __cplusplus
}
#endif


#endif
