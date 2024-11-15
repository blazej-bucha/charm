/* This header file is not a part of API. */


#ifndef __CMP_VALS_H__
#define __CMP_VALS_H__


#include <config.h>
#include "../src/prec.h"


#ifdef __cplusplus
extern "C"
{
#endif


extern long int cmp_vals_real(REAL,
                              REAL,
                              REAL);


extern long int cmp_vals_ulong(unsigned long,
                               unsigned long);


#ifdef __cplusplus
}
#endif


#endif
